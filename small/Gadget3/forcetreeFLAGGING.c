#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

#ifdef SUBFIND
#include "subfind.h"
#endif

#ifdef NUM_THREADS
#include <pthread.h>
#endif


/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */

/*! auxialiary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float shortrange_table[NTAB], shortrange_table_potential[NTAB];
#ifdef DISTORTIONTENSORPS
static float shortrange_table_tidal[NTAB];
#endif
/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;

static int tree_allocated_flag = 0;


#ifdef NUM_THREADS
extern pthread_mutex_t mutex_nexport, mutex_partnodedrift, mutex_workcount;

#define LOCK_NEXPORT         pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT       pthread_mutex_unlock(&mutex_nexport);
#define LOCK_PARTNODEDRIFT   pthread_mutex_lock(&mutex_partnodedrift);
#define UNLOCK_PARTNODEDRIFT pthread_mutex_unlock(&mutex_partnodedrift);

/*! The cost computation for the tree-gravity (required for the domain
decomposition) is not exactly thread-safe if THREAD_SAFE_COSTS is not defined. 
However using locks for an exactly thread-safe cost computiation results in a
significant (~25%) performance penalty in the tree-walk while having only an 
extremely small effect on the obtained costs. The domain decomposition should
thus not be significantly changed if THREAD_SAFE_COSTS is not used.*/
#ifdef THREAD_SAFE_COSTS
#define LOCK_WORKCOUNT       pthread_mutex_lock(&mutex_workcount);
#define UNLOCK_WORKCOUNT     pthread_mutex_unlock(&mutex_workcount);
#else
#define LOCK_WORKCOUNT
#define UNLOCK_WORKCOUNT
#endif

#else

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT
#define LOCK_WORKCOUNT
#define UNLOCK_WORKCOUNT

#endif



#ifdef PERIODIC
/*! Size of 3D lock-up table for Ewald correction force */
#define EN  64
/*! 3D lock-up table for Ewald correction to force and potential. Only one
 *  octant is stored, the rest constructed by using the symmetry
 */
static MyFloat fcorrx[EN + 1][EN + 1][EN + 1];
static MyFloat fcorry[EN + 1][EN + 1][EN + 1];
static MyFloat fcorrz[EN + 1][EN + 1][EN + 1];
static MyFloat potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
#endif



/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart, struct unbind_data *mp)
{
	int flag;

	do
	{
		Numnodestree = force_treebuild_single(npart, mp);

		MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		if(flag == -1)
		{
			force_treefree();

			if(ThisTask == 0)
				printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);

			All.TreeAllocFactor *= 1.15;

			if(ThisTask == 0)
				printf("new value=%g\n", All.TreeAllocFactor);

			force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
		}
	}
	while(flag == -1);

	force_flag_localnodes();

	force_exchange_pseudodata();

	force_treeupdate_pseudos(All.MaxPart);

	TimeOfLastTreeConstruction = All.Time;

	return Numnodestree;
}



/*! Constructs the gravitational oct-tree.  
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart, struct unbind_data *mp)
{
	int i, j, k, subnode = 0, shift, parent, numnodes, rep;
	int nfree, th, nn, no;
	struct NODE *nfreep;
	MyFloat lenhalf;
	peanokey key, morton, th_key, *morton_list;


	/* create an empty root node  */
	nfree = All.MaxPart;		/* index of first free node */
	nfreep = &Nodes[nfree];	/* select first node */

	nfreep->len = DomainLen;
	for(j = 0; j < 3; j++)
		nfreep->center[j] = DomainCenter[j];
	for(j = 0; j < 8; j++)
		nfreep->u.suns[j] = -1;


	numnodes = 1;
	nfreep++;
	nfree++;

	/* create a set of empty nodes corresponding to the top-level domain
	 * grid. We need to generate these nodes first to make sure that we have a
	 * complete top-level tree which allows the easy insertion of the
	 * pseudo-particles at the right place
	 */

	force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);

	/* if a high-resolution region in a global tree is used, we need to generate
	 * an additional set empty nodes to make sure that we have a complete
	 * top-level tree for the high-resolution inset
	 */

	nfreep = &Nodes[nfree];
	parent = -1;			/* note: will not be used below before it is changed */

	morton_list = (peanokey *) mymalloc("morton_list", NumPart * sizeof(peanokey));

	/* now we insert all particles */
	for(k = 0; k < npart; k++)
	{
		if(mp)
			i = mp[k].index;
		else
			i = k;

#ifdef NEUTRINOS
		if(P[i].Type == 2)
			continue;
#endif

		rep = 0;

		key = peano_and_morton_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
				(int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
				(int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION,
				&morton);
		morton_list[i] = morton;

		shift = 3 * (BITS_PER_DIMENSION - 1);

		no = 0;
		while(TopNodes[no].Daughter >= 0)
		{
			no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);
			shift -= 3;
			rep++;
		}

		no = TopNodes[no].Leaf;
		th = DomainNodeIndex[no];

		while(1)
		{
			if(th >= All.MaxPart)	/* we are dealing with an internal node */
			{
				if(shift >= 0)
				{
					subnode = ((morton >> shift) & 7);
				}
				else
				{
					subnode = 0;
					if(P[i].Pos[0] > Nodes[th].center[0])
						subnode += 1;
					if(P[i].Pos[1] > Nodes[th].center[1])
						subnode += 2;
					if(P[i].Pos[2] > Nodes[th].center[2])
						subnode += 4;
				}

#ifndef NOTREERND
				if(Nodes[th].len < 1.0e-3 * All.ForceSoftening[P[i].Type])
				{
					/* seems like we're dealing with particles at identical (or extremely close)
					 * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
					 * of tree are still correct, but this will only happen well below gravitational softening
					 * length-scale anyway.
					 */
					subnode = (int) (8.0 * get_random_number((P[i].ID + rep) % (RNDTABLE + (rep & 3))));

					if(subnode >= 8)
						subnode = 7;
				}
#endif

				nn = Nodes[th].u.suns[subnode];

				shift -= 3;

				if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
				{
					parent = th;
					th = nn;
					rep++;
				}
				else
				{
					/* here we have found an empty slot where we can attach
					 * the new particle as a leaf.
					 */
					Nodes[th].u.suns[subnode] = i;
					break;	/* done for this particle */
				}
			}
			else
			{
				/* We try to insert into a leaf with a single particle.  Need
				 * to generate a new internal node at this point.
				 */
				Nodes[parent].u.suns[subnode] = nfree;

				nfreep->len = 0.5 * Nodes[parent].len;
				lenhalf = 0.25 * Nodes[parent].len;

				if(subnode & 1)
					nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
				else
					nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

				if(subnode & 2)
					nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
				else
					nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

				if(subnode & 4)
					nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
				else
					nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

				nfreep->u.suns[0] = -1;
				nfreep->u.suns[1] = -1;
				nfreep->u.suns[2] = -1;
				nfreep->u.suns[3] = -1;
				nfreep->u.suns[4] = -1;
				nfreep->u.suns[5] = -1;
				nfreep->u.suns[6] = -1;
				nfreep->u.suns[7] = -1;

				if(shift >= 0)
				{
					th_key = morton_list[th];
					subnode = ((th_key >> shift) & 7);
				}
				else
				{
					subnode = 0;
					if(P[th].Pos[0] > nfreep->center[0])
						subnode += 1;
					if(P[th].Pos[1] > nfreep->center[1])
						subnode += 2;
					if(P[th].Pos[2] > nfreep->center[2])
						subnode += 4;
				}

#ifndef NOTREERND
				if(nfreep->len < 1.0e-3 * All.ForceSoftening[P[th].Type])
				{
					/* seems like we're dealing with particles at identical (or extremely close)
					 * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
					 * of tree are still correct, but this will only happen well below gravitational softening
					 * length-scale anyway.
					 */
					subnode = (int) (8.0 * get_random_number((P[th].ID + rep) % (RNDTABLE + (rep & 3))));

					if(subnode >= 8)
						subnode = 7;
				}
#endif
				nfreep->u.suns[subnode] = th;

				th = nfree;	/* resume trying to insert the new particle at
				 * the newly created internal node
				 */

				numnodes++;
				nfree++;
				nfreep++;

				if((numnodes) >= MaxNodes)
				{
					printf("task %d: maximum number %d of tree-nodes reached for particle %d.\n", ThisTask,
							MaxNodes, i);

					if(All.TreeAllocFactor > 5.0)
					{
						printf
						("task %d: looks like a serious problem for particle %d, stopping with particle dump.\n",
								ThisTask, i);
						dump_particles();
						endrun(1);
					}
					else
					{
						myfree(morton_list);
						return -1;
					}
				}
			}
		}
	}

	myfree(morton_list);


	/* insert the pseudo particles that represent the mass distribution of other domains */
	force_insert_pseudo_particles();


	/* now compute the multipole moments recursively */
	last = -1;

	force_update_node_recursive(All.MaxPart, -1, -1);

	if(last >= All.MaxPart)
	{
		if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
			Nextnode[last - MaxNodes] = -1;
		else
			Nodes[last].u.d.nextnode = -1;
	}
	else
		Nextnode[last] = -1;

	return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
		int *nextfree)
{
	int i, j, k, n, sub, count;
	MyFloat lenhalf;

	if(TopNodes[topnode].Daughter >= 0)
	{
		for(i = 0; i < 2; i++)
			for(j = 0; j < 2; j++)
				for(k = 0; k < 2; k++)
				{
					sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

					count = i + 2 * j + 4 * k;

					Nodes[no].u.suns[count] = *nextfree;

					lenhalf = 0.25 * Nodes[no].len;
					Nodes[*nextfree].len = 0.5 * Nodes[no].len;
					Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
					Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
					Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

					for(n = 0; n < 8; n++)
						Nodes[*nextfree].u.suns[n] = -1;

					if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
						DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

					*nextfree = *nextfree + 1;
					*nodecount = *nodecount + 1;

					if((*nodecount) >= MaxNodes || (*nodecount) >= MaxTopNodes)
					{
						printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
								"MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
								ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
						printf("in create empty nodes\n");
						dump_particles();
						endrun(11);
					}

					force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
							bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
				}
	}
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
	int i, index;

	for(i = 0; i < NTopleaves; i++)
	{
		index = DomainNodeIndex[i];

		if(DomainTask[i] != ThisTask)
			Nodes[index].u.suns[0] = All.MaxPart + MaxNodes + i;
	}
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 * 
 *  If UNEQUALSOFTENINGS is set, bits 2-4 give the particle type with
 *  the maximum softening among the particles in the node, and bit 5
 *  flags whether the node contains any particles with lower softening
 *  than that.  
 */
void force_update_node_recursive(int no, int sib, int father)
{
	int j, jj, k, p, pp, nextsib, suns[8], count_particles, multiple_flag;
	MyFloat hmax, vmax, v;
#ifdef KD_HMAX_ESTIMATE
	integertime ti_step, time0, time1;
	double dt_drift_hmax;
#else
	MyFloat divVmax;
#endif
	MyFloat s[3], vs[3], mass;
	struct particle_data *pa;

#ifdef SCALARFIELD
	MyFloat s_dm[3], vs_dm[3], mass_dm;
#endif
#ifdef RADTRANSFER
	MyFloat stellar_mass;
	MyFloat stellar_s[3];
#endif
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	int maxsofttype, current_maxsofttype, diffsoftflag;
#else
	MyFloat maxsoft;
#endif
#endif

#ifdef ADAPTGRAVSOFT
	MyFloat maxsoft;
#endif

	if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
	{
		for(j = 0; j < 8; j++)
#ifdef GRAVITY_CENTROID
			suns[j] = Extnodes[no].suns[j] = Nodes[no].u.suns[j];
#else
		suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will overwrite one element (union!) */
#endif
		if(last >= 0)
		{
			if(last >= All.MaxPart)
			{
				if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
					Nextnode[last - MaxNodes] = no;
				else
					Nodes[last].u.d.nextnode = no;
			}
			else
				Nextnode[last] = no;
		}

		last = no;

#ifdef RADTRANSFER
		stellar_mass = 0;
		stellar_s[0] = 0;
		stellar_s[1] = 0;
		stellar_s[2] = 0;
#endif
#ifdef SCALARFIELD
		mass_dm = 0;
		s_dm[0] = vs_dm[0] = 0;
		s_dm[1] = vs_dm[1] = 0;
		s_dm[2] = vs_dm[2] = 0;
#endif
		mass = 0;
		s[0] = 0;
		s[1] = 0;
		s[2] = 0;
		vs[0] = 0;
		vs[1] = 0;
		vs[2] = 0;
		hmax = 0;
		vmax = 0;
#ifndef KD_HMAX_ESTIMATE
		divVmax = 0;
#endif
		count_particles = 0;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		maxsofttype = 7;
		diffsoftflag = 0;
#else
		maxsoft = 0;
#endif
#endif

#ifdef ADAPTGRAVSOFT
		maxsoft = 0;
#endif

		for(j = 0; j < 8; j++)
		{
			if((p = suns[j]) >= 0)
			{
				/* check if we have a sibling on the same level */
				for(jj = j + 1; jj < 8; jj++)
					if((pp = suns[jj]) >= 0)
						break;

				if(jj < 8)	/* yes, we do */
					nextsib = pp;
				else
					nextsib = sib;

				force_update_node_recursive(p, nextsib, no);

				if(p >= All.MaxPart)	/* an internal node or pseudo particle */
				{
					if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
					{
						/* nothing to be done here because the mass of the
						 * pseudo-particle is still zero. This will be changed
						 * later.
						 */
					}
					else
					{
						mass += (Nodes[p].u.d.mass);
						s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
						s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
						s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
						vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
						vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
						vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);
#ifdef RADTRANSFER
						stellar_s[0] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[0]);
						stellar_s[1] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[1]);
						stellar_s[2] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[2]);
						stellar_mass += (Nodes[p].stellar_mass);
#endif
#ifdef SCALARFIELD
						mass_dm += (Nodes[p].mass_dm);
						s_dm[0] += (Nodes[p].mass_dm * Nodes[p].s_dm[0]);
						s_dm[1] += (Nodes[p].mass_dm * Nodes[p].s_dm[1]);
						s_dm[2] += (Nodes[p].mass_dm * Nodes[p].s_dm[2]);
						vs_dm[0] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[0]);
						vs_dm[1] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[1]);
						vs_dm[2] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[2]);
#endif
						if(Nodes[p].u.d.mass > 0)
						{
							if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
								count_particles += 2;
							else
								count_particles++;
						}

						if(Extnodes[p].hmax > hmax)
							hmax = Extnodes[p].hmax;

						if(Extnodes[p].vmax > vmax)
							vmax = Extnodes[p].vmax;
#ifndef KD_HMAX_ESTIMATE
						if(Extnodes[p].divVmax > divVmax)
							divVmax = Extnodes[p].divVmax;
#endif

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
						diffsoftflag |= maskout_different_softening_flag(Nodes[p].u.d.bitflags);

						if(maxsofttype == 7)
							maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
						else
						{
							current_maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
							if(current_maxsofttype != 7)
							{
								if(All.ForceSoftening[current_maxsofttype] > All.ForceSoftening[maxsofttype])
								{
									maxsofttype = current_maxsofttype;
									diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
								}
								else
								{
									if(All.ForceSoftening[current_maxsofttype] <
											All.ForceSoftening[maxsofttype])
										diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
								}
							}
						}
#else
						if(Nodes[p].maxsoft > maxsoft)
							maxsoft = Nodes[p].maxsoft;
#endif
#endif

#ifdef ADAPTGRAVSOFT
						/* update of the maximum gravitational softening in the node */
						if(Nodes[p].maxsoft > maxsoft)
							maxsoft = Nodes[p].maxsoft;
#endif

					}
				}
				else		/* a particle */
				{
					count_particles++;

					pa = &P[p];

					mass += (pa->Mass);
#ifdef GRAVITY_CENTROID
					if(P[p].Type == 0)
					{
						s[0] += (pa->Mass * SphP[p].Center[0]);
						s[1] += (pa->Mass * SphP[p].Center[1]);
						s[2] += (pa->Mass * SphP[p].Center[2]);
					}
					else
					{
						s[0] += (pa->Mass * pa->Pos[0]);
						s[1] += (pa->Mass * pa->Pos[1]);
						s[2] += (pa->Mass * pa->Pos[2]);
					}
#else
	s[0] += (pa->Mass * pa->Pos[0]);
	s[1] += (pa->Mass * pa->Pos[1]);
	s[2] += (pa->Mass * pa->Pos[2]);
#endif
	vs[0] += (pa->Mass * pa->Vel[0]);
	vs[1] += (pa->Mass * pa->Vel[1]);
	vs[2] += (pa->Mass * pa->Vel[2]);
#ifdef RADTRANSFER
#ifdef EDDINGTON_TENSOR_STARS
	if(pa->Type == 4)
	{
		stellar_s[0] += (pa->Mass * pa->Pos[0]);
		stellar_s[1] += (pa->Mass * pa->Pos[1]);
		stellar_s[2] += (pa->Mass * pa->Pos[2]);
		stellar_mass += (pa->Mass);
	}
#endif
#ifdef EDDINGTON_TENSOR_GAS
	if(pa->Type == 0)
	{
		stellar_s[0] += (pa->Mass * pa->Pos[0]);
		stellar_s[1] += (pa->Mass * pa->Pos[1]);
		stellar_s[2] += (pa->Mass * pa->Pos[2]);
		stellar_mass += (pa->Mass);
	}
#endif

#if defined(SFR) && defined(EDDINGTON_TENSOR_SFR)
	if(pa->Type == 0)
	{
		double a;

		if(All.ComovingIntegrationOn == 1)
			a = All.Time;
		else
			a = 1.0;

		if((SphP[p].d.Density / (a * a * a)) >= All.PhysDensThresh)
		{
			stellar_s[0] += (pa->Mass * pa->Pos[0]);
			stellar_s[1] += (pa->Mass * pa->Pos[1]);
			stellar_s[2] += (pa->Mass * pa->Pos[2]);
			stellar_mass += (pa->Mass);
		}
	}
#endif
#if defined(BLACK_HOLES) && defined(EDDINGTON_TENSOR_BH)
	if(pa->Type == 5)
	{
		stellar_s[0] += (pa->Mass * pa->Pos[0]);
		stellar_s[1] += (pa->Mass * pa->Pos[1]);
		stellar_s[2] += (pa->Mass * pa->Pos[2]);
		stellar_mass += (pa->Mass);
	}
#endif
#endif

#ifdef SCALARFIELD
	if(pa->Type != 0)
	{
		mass_dm += (pa->Mass);
		s_dm[0] += (pa->Mass * pa->Pos[0]);
		s_dm[1] += (pa->Mass * pa->Pos[1]);
		s_dm[2] += (pa->Mass * pa->Pos[2]);
		vs_dm[0] += (pa->Mass * pa->Vel[0]);
		vs_dm[1] += (pa->Mass * pa->Vel[1]);
		vs_dm[2] += (pa->Mass * pa->Vel[2]);
	}
#endif
if(pa->Type == 0)
{
#ifdef KD_HMAX_ESTIMATE
	ti_step = P[p].TimeBin ? (((integertime) 1) << P[p].TimeBin) : 0;
	time0 = All.Ti_Current;
	time1 = P[p].Ti_begstep + ti_step;
	dt_drift_hmax = get_drift_factor(time0, time1);

	hmax = DMAX(hmax, PPP[p].Hsml);
	hmax = DMAX(hmax,
			(PPP[p].Hsml * exp(0.333333333333 * SphP[p].v.DivVel * dt_drift_hmax)));
#else
	if(PPP[p].Hsml > hmax)
		hmax = PPP[p].Hsml;

	if(SphP[p].v.DivVel > divVmax)
		divVmax = SphP[p].v.DivVel;
#endif
}

for(k = 0; k < 3; k++)
	if((v = fabs(pa->Vel[k])) > vmax)
		vmax = v;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
if(maxsofttype == 7)
	maxsofttype = pa->Type;
else
{
	if(All.ForceSoftening[pa->Type] > All.ForceSoftening[maxsofttype])
	{
		maxsofttype = pa->Type;
		diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
	}
	else
	{
		if(All.ForceSoftening[pa->Type] < All.ForceSoftening[maxsofttype])
			diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
	}
}
#else
	if(pa->Type == 0)
	{
		if(PPP[p].Hsml > maxsoft)
			maxsoft = PPP[p].Hsml;
	}
	else
	{
		if(All.ForceSoftening[pa->Type] > maxsoft)
			maxsoft = All.ForceSoftening[pa->Type];
	}
#endif
#endif

#ifdef ADAPTGRAVSOFT
	/* update of the maximum gravitational softening  */
	if(pa->AGS_Hsml > maxsoft)
		maxsoft = pa->AGS_Hsml;
#endif
				}
			}
		}


		if(mass)
		{
			s[0] /= mass;
			s[1] /= mass;
			s[2] /= mass;
			vs[0] /= mass;
			vs[1] /= mass;
			vs[2] /= mass;
		}
		else
		{
			s[0] = Nodes[no].center[0];
			s[1] = Nodes[no].center[1];
			s[2] = Nodes[no].center[2];
			vs[0] = 0;
			vs[1] = 0;
			vs[2] = 0;
		}

#ifdef RADTRANSFER
		if(stellar_mass)
		{
			stellar_s[0] /= stellar_mass;
			stellar_s[1] /= stellar_mass;
			stellar_s[2] /= stellar_mass;
		}
		else
		{
			stellar_s[0] = Nodes[no].center[0];
			stellar_s[1] = Nodes[no].center[1];
			stellar_s[2] = Nodes[no].center[2];
		}
#endif
#ifdef SCALARFIELD
		if(mass_dm)
		{
			s_dm[0] /= mass_dm;
			s_dm[1] /= mass_dm;
			s_dm[2] /= mass_dm;
			vs_dm[0] /= mass_dm;
			vs_dm[1] /= mass_dm;
			vs_dm[2] /= mass_dm;
		}
		else
		{
			s_dm[0] = Nodes[no].center[0];
			s_dm[1] = Nodes[no].center[1];
			s_dm[2] = Nodes[no].center[2];
			vs_dm[0] = 0;
			vs_dm[1] = 0;
			vs_dm[2] = 0;
		}
#endif


		Nodes[no].Ti_current = All.Ti_Current;
		Nodes[no].u.d.mass = mass;
		Nodes[no].u.d.s[0] = s[0];
		Nodes[no].u.d.s[1] = s[1];
		Nodes[no].u.d.s[2] = s[2];
		Nodes[no].GravCost = 0;
#ifdef RADTRANSFER
		Nodes[no].stellar_s[0] = stellar_s[0];
		Nodes[no].stellar_s[1] = stellar_s[1];
		Nodes[no].stellar_s[2] = stellar_s[2];
		Nodes[no].stellar_mass = stellar_mass;
#endif
#ifdef SCALARFIELD
		Nodes[no].s_dm[0] = s_dm[0];
		Nodes[no].s_dm[1] = s_dm[1];
		Nodes[no].s_dm[2] = s_dm[2];
		Nodes[no].mass_dm = mass_dm;
		Extnodes[no].vs_dm[0] = vs_dm[0];
		Extnodes[no].vs_dm[1] = vs_dm[1];
		Extnodes[no].vs_dm[2] = vs_dm[2];
		Extnodes[no].dp_dm[0] = 0;
		Extnodes[no].dp_dm[1] = 0;
		Extnodes[no].dp_dm[2] = 0;
#endif

Extnodes[no].Ti_lastkicked = All.Ti_Current;
Extnodes[no].Flag = GlobFlag;
Extnodes[no].vs[0] = vs[0];
Extnodes[no].vs[1] = vs[1];
Extnodes[no].vs[2] = vs[2];
Extnodes[no].hmax = hmax;
Extnodes[no].vmax = vmax;
#ifndef KD_HMAX_ESTIMATE
Extnodes[no].divVmax = divVmax;
#endif
Extnodes[no].dp[0] = 0;
Extnodes[no].dp[1] = 0;
Extnodes[no].dp[2] = 0;

#ifdef FLTROUNDOFFREDUCTION
Extnodes[no].s_base[0] = Nodes[no].u.d.s[0];
Extnodes[no].s_base[1] = Nodes[no].u.d.s[1];
Extnodes[no].s_base[2] = Nodes[no].u.d.s[2];
Extnodes[no].len_base = Nodes[no].len;
#ifdef SCALARFIELD
Extnodes[no].s_dm_base[0] = Nodes[no].s_dm[0];
Extnodes[no].s_dm_base[1] = Nodes[no].s_dm[1];
Extnodes[no].s_dm_base[2] = Nodes[no].s_dm[2];
#endif
#endif
if(count_particles > 1)	/* this flags that the node represents more than one particle */
	multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
else
	multiple_flag = 0;

Nodes[no].u.d.bitflags = multiple_flag;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
Nodes[no].u.d.bitflags |= diffsoftflag + (maxsofttype << BITFLAG_MAX_SOFTENING_TYPE);
#else
Nodes[no].maxsoft = maxsoft;
#endif
#endif
#ifdef ADAPTGRAVSOFT
Nodes[no].maxsoft = maxsoft;
#endif
Nodes[no].u.d.sibling = sib;
Nodes[no].u.d.father = father;
	}
	else				/* single particle or pseudo particle */
	{
		if(last >= 0)
		{
			if(last >= All.MaxPart)
			{
				if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
					Nextnode[last - MaxNodes] = no;
				else
					Nodes[last].u.d.nextnode = no;
			}
			else
				Nextnode[last] = no;
		}

		last = no;

		if(no < All.MaxPart)	/* only set it for single particles */
			Father[no] = father;
	}
}




/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
	int i, no, m, ta, recvTask;
	int *recvcounts, *recvoffset;
	struct DomainNODE
	{
		MyFloat s[3];
		MyFloat vs[3];
		MyFloat mass;
		MyFloat hmax;
		MyFloat vmax;
#ifndef KD_HMAX_ESTIMATE
		MyFloat divVmax;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTGRAVSOFT)
		MyFloat maxsoft;
#endif
#ifdef RADTRANSFER
		MyFloat stellar_mass;
		MyFloat stellar_s[3];
#endif
#ifdef SCALARFIELD
		MyFloat s_dm[3];
		MyFloat vs_dm[3];
		MyFloat mass_dm;
#endif
		unsigned int bitflags;
#ifdef PAD_STRUCTURES
#ifndef DOUBLEPRECISION
		int pad[5];
#else
#if (DOUBLEPRECISION+0) == 2
		/* mixed precision */
		int pad[5];
#else
		int pad[3];
#endif
#endif				/* DOUBLEPRECISION  */
#endif
	}
	*DomainMoment;


	DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

	for(m = 0; m < MULTIPLEDOMAINS; m++)
		for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
				i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
		{
			no = DomainNodeIndex[i];

			/* read out the multipole moments from the local base cells */
			DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
			DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
			DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
			DomainMoment[i].vs[0] = Extnodes[no].vs[0];
			DomainMoment[i].vs[1] = Extnodes[no].vs[1];
			DomainMoment[i].vs[2] = Extnodes[no].vs[2];
			DomainMoment[i].mass = Nodes[no].u.d.mass;
			DomainMoment[i].hmax = Extnodes[no].hmax;
			DomainMoment[i].vmax = Extnodes[no].vmax;
#ifndef KD_HMAX_ESTIMATE
			DomainMoment[i].divVmax = Extnodes[no].divVmax;
#endif
			DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTGRAVSOFT)
			DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
#ifdef RADTRANSFER
			DomainMoment[i].stellar_s[0] = Nodes[no].stellar_s[0];
			DomainMoment[i].stellar_s[1] = Nodes[no].stellar_s[1];
			DomainMoment[i].stellar_s[2] = Nodes[no].stellar_s[2];
			DomainMoment[i].stellar_mass = Nodes[no].stellar_mass;
#endif
#ifdef SCALARFIELD
			DomainMoment[i].s_dm[0] = Nodes[no].s_dm[0];
			DomainMoment[i].s_dm[1] = Nodes[no].s_dm[1];
			DomainMoment[i].s_dm[2] = Nodes[no].s_dm[2];
			DomainMoment[i].mass_dm = Nodes[no].mass_dm;
			DomainMoment[i].vs_dm[0] = Extnodes[no].vs_dm[0];
			DomainMoment[i].vs_dm[1] = Extnodes[no].vs_dm[1];
			DomainMoment[i].vs_dm[2] = Extnodes[no].vs_dm[2];
#endif
		}

	/* share the pseudo-particle data accross CPUs */

	recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
	recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);

	for(m = 0; m < MULTIPLEDOMAINS; m++)
	{
		for(recvTask = 0; recvTask < NTask; recvTask++)
		{
			recvcounts[recvTask] =
					(DomainEndList[recvTask * MULTIPLEDOMAINS + m] - DomainStartList[recvTask * MULTIPLEDOMAINS + m] +
							1) * sizeof(struct DomainNODE);
			recvoffset[recvTask] = DomainStartList[recvTask * MULTIPLEDOMAINS + m] * sizeof(struct DomainNODE);
		}

		MPI_Allgatherv(&DomainMoment[DomainStartList[ThisTask * MULTIPLEDOMAINS + m]], recvcounts[ThisTask],
				MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
	}

	myfree(recvoffset);
	myfree(recvcounts);


	for(ta = 0; ta < NTask; ta++)
		if(ta != ThisTask)
			for(m = 0; m < MULTIPLEDOMAINS; m++)
				for(i = DomainStartList[ta * MULTIPLEDOMAINS + m]; i <= DomainEndList[ta * MULTIPLEDOMAINS + m]; i++)
				{
					no = DomainNodeIndex[i];

					Nodes[no].u.d.s[0] = DomainMoment[i].s[0];
					Nodes[no].u.d.s[1] = DomainMoment[i].s[1];
					Nodes[no].u.d.s[2] = DomainMoment[i].s[2];
#ifdef FLTROUNDOFFREDUCTION
					Extnodes[no].s_base[0] = DomainMoment[i].s[0];
					Extnodes[no].s_base[1] = DomainMoment[i].s[1];
					Extnodes[no].s_base[2] = DomainMoment[i].s[2];
#endif
					Extnodes[no].vs[0] = DomainMoment[i].vs[0];
					Extnodes[no].vs[1] = DomainMoment[i].vs[1];
					Extnodes[no].vs[2] = DomainMoment[i].vs[2];
					Nodes[no].u.d.mass = DomainMoment[i].mass;
					Extnodes[no].hmax = DomainMoment[i].hmax;
					Extnodes[no].vmax = DomainMoment[i].vmax;
#ifndef KD_HMAX_ESTIMATE
					Extnodes[no].divVmax = DomainMoment[i].divVmax;
#endif
					Nodes[no].u.d.bitflags =
							(Nodes[no].u.d.bitflags & (~BITFLAG_MASK)) | (DomainMoment[i].bitflags & BITFLAG_MASK);
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTGRAVSOFT)
					Nodes[no].maxsoft = DomainMoment[i].maxsoft;
#endif
#ifdef RADTRANSFER
					Nodes[no].stellar_s[0] = DomainMoment[i].stellar_s[0];
					Nodes[no].stellar_s[1] = DomainMoment[i].stellar_s[1];
					Nodes[no].stellar_s[2] = DomainMoment[i].stellar_s[2];
					Nodes[no].stellar_mass = DomainMoment[i].stellar_mass;
#endif
#ifdef SCALARFIELD
					Nodes[no].s_dm[0] = DomainMoment[i].s_dm[0];
					Nodes[no].s_dm[1] = DomainMoment[i].s_dm[1];
					Nodes[no].s_dm[2] = DomainMoment[i].s_dm[2];
					Nodes[no].mass_dm = DomainMoment[i].mass_dm;
					Extnodes[no].vs_dm[0] = DomainMoment[i].vs_dm[0];
					Extnodes[no].vs_dm[1] = DomainMoment[i].vs_dm[1];
					Extnodes[no].vs_dm[2] = DomainMoment[i].vs_dm[2];
#endif
				}

	myfree(DomainMoment);
}



#ifdef GRAVITY_CENTROID
void force_update_node_center_of_mass_recursive(int no, int sib, int father)
{
	int j, jj, p, pp, nextsib, suns[8], count_particles;

	//   int k, multiple_flag;
	//   MyFloat hmax, vmax, v, divVmax;
	MyFloat s[3], mass;

	//   MyFloat vs[3];
	struct particle_data *pa;


	if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
	{
		for(j = 0; j < 8; j++)
			suns[j] = Extnodes[no].suns[j];	/* this "backup" is necessary because the nextnode entry will */

		mass = 0;
		s[0] = 0;
		s[1] = 0;
		s[2] = 0;

		for(j = 0; j < 8; j++)
		{
			if((p = suns[j]) >= 0)
			{
				/* check if we have a sibling on the same level */
				for(jj = j + 1; jj < 8; jj++)
					if((pp = suns[jj]) >= 0)
						break;

				if(jj < 8)	/* yes, we do */
					nextsib = pp;
				else
					nextsib = sib;

				force_update_node_center_of_mass_recursive(p, nextsib, no);

				if(p >= All.MaxPart)	/* an internal node or pseudo particle */
				{
					if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
					{
						/* nothing to be done here because the mass of the
						 * pseudo-particle is still zero. This will be changed
						 * later.
						 */
					}
					else
					{
						mass += (Nodes[p].u.d.mass);
						s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
						s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
						s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);


					}
				}
				else		/* a particle */
				{
					count_particles++;

					pa = &P[p];

					mass += (pa->Mass);

					if(P[p].Type == 0)
					{
						s[0] += (pa->Mass * SphP[p].Center[0]);
						s[1] += (pa->Mass * SphP[p].Center[1]);
						s[2] += (pa->Mass * SphP[p].Center[2]);
					}
					else
					{
						s[0] += (pa->Mass * pa->Pos[0]);
						s[1] += (pa->Mass * pa->Pos[1]);
						s[2] += (pa->Mass * pa->Pos[2]);
					}
				}
			}
		}


		if(mass)
		{
			s[0] /= mass;
			s[1] /= mass;
			s[2] /= mass;
		}
		else
		{
			s[0] = Nodes[no].center[0];
			s[1] = Nodes[no].center[1];
			s[2] = Nodes[no].center[2];
		}

		Nodes[no].u.d.s[0] = s[0];
		Nodes[no].u.d.s[1] = s[1];
		Nodes[no].u.d.s[2] = s[2];
	}
}
#endif


/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
	int j, p, count_particles, multiple_flag;
	MyFloat hmax, vmax;
#ifndef KD_HMAX_ESTIMATE
	MyFloat divVmax;
#endif
	MyFloat s[3], vs[3], mass;

#ifdef RADTRANSFER
	MyFloat stellar_mass;
	MyFloat stellar_s[3];

#endif
#ifdef SCALARFIELD
	MyFloat s_dm[3], vs_dm[3], mass_dm;
#endif

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	int maxsofttype, diffsoftflag, current_maxsofttype;
#else
	MyFloat maxsoft;
#endif
#endif

#ifdef ADAPTGRAVSOFT
	MyFloat maxsoft;
#endif

#ifdef RADTRANSFER
	stellar_mass = 0;
	stellar_s[0] = 0;
	stellar_s[1] = 0;
	stellar_s[2] = 0;
#endif
#ifdef SCALARFIELD
	mass_dm = 0;
	s_dm[0] = vs_dm[0] = 0;
	s_dm[1] = vs_dm[1] = 0;
	s_dm[2] = vs_dm[2] = 0;
#endif
	mass = 0;
	s[0] = 0;
	s[1] = 0;
	s[2] = 0;
	vs[0] = 0;
	vs[1] = 0;
	vs[2] = 0;
	hmax = 0;
	vmax = 0;
#ifndef KD_HMAX_ESTIMATE
	divVmax = 0;
#endif
	count_particles = 0;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	maxsofttype = 7;
	diffsoftflag = 0;
#else
	maxsoft = 0;
#endif
#endif

#ifdef ADAPTGRAVSOFT
	maxsoft = 0;
#endif

	p = Nodes[no].u.d.nextnode;

	for(j = 0; j < 8; j++)	/* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
	{
		if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)	/* internal node */
		{
			if(Nodes[p].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
				force_treeupdate_pseudos(p);

			mass += (Nodes[p].u.d.mass);
			s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
			s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
			s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
#ifdef RADTRANSFER
			stellar_mass += (Nodes[p].stellar_mass);
			stellar_s[0] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[0]);
			stellar_s[1] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[1]);
			stellar_s[2] += (Nodes[p].stellar_mass * Nodes[p].stellar_s[2]);
#endif
#ifdef SCALARFIELD
			mass_dm += (Nodes[p].mass_dm);
			s_dm[0] += (Nodes[p].mass_dm * Nodes[p].s_dm[0]);
			s_dm[1] += (Nodes[p].mass_dm * Nodes[p].s_dm[1]);
			s_dm[2] += (Nodes[p].mass_dm * Nodes[p].s_dm[2]);
			vs_dm[0] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[0]);
			vs_dm[1] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[1]);
			vs_dm[2] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[2]);
#endif
			vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
			vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
			vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);

			if(Extnodes[p].hmax > hmax)
				hmax = Extnodes[p].hmax;
			if(Extnodes[p].vmax > vmax)
				vmax = Extnodes[p].vmax;
#ifndef KD_HMAX_ESTIMATE
			if(Extnodes[p].divVmax > divVmax)
				divVmax = Extnodes[p].divVmax;
#endif

			if(Nodes[p].u.d.mass > 0)
			{
				if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
					count_particles += 2;
				else
					count_particles++;
			}

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
			diffsoftflag |= maskout_different_softening_flag(Nodes[p].u.d.bitflags);

			if(maxsofttype == 7)
				maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
			else
			{
				current_maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
				if(current_maxsofttype != 7)
				{
					if(All.ForceSoftening[current_maxsofttype] > All.ForceSoftening[maxsofttype])
					{
						maxsofttype = current_maxsofttype;
						diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
					}
					else
					{
						if(All.ForceSoftening[current_maxsofttype] < All.ForceSoftening[maxsofttype])
							diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
					}
				}
			}
#else
	if(Nodes[p].maxsoft > maxsoft)
		maxsoft = Nodes[p].maxsoft;
#endif
#endif

#ifdef ADAPTGRAVSOFT
	if(Nodes[p].maxsoft > maxsoft)
		maxsoft = Nodes[p].maxsoft;
#endif
		}
		else
			endrun(6767);		/* may not happen */

		p = Nodes[p].u.d.sibling;
	}

	if(mass)
	{
		s[0] /= mass;
		s[1] /= mass;
		s[2] /= mass;
		vs[0] /= mass;
		vs[1] /= mass;
		vs[2] /= mass;
	}
	else
	{
		s[0] = Nodes[no].center[0];
		s[1] = Nodes[no].center[1];
		s[2] = Nodes[no].center[2];
		vs[0] = 0;
		vs[1] = 0;
		vs[2] = 0;
	}

#ifdef RADTRANSFER
	if(stellar_mass)
	{
		stellar_s[0] /= stellar_mass;
		stellar_s[1] /= stellar_mass;
		stellar_s[2] /= stellar_mass;
	}
	else
	{
		stellar_s[0] = Nodes[no].center[0];
		stellar_s[1] = Nodes[no].center[1];
		stellar_s[2] = Nodes[no].center[2];
	}
#endif
#ifdef SCALARFIELD
if(mass_dm)
{
	s_dm[0] /= mass_dm;
	s_dm[1] /= mass_dm;
	s_dm[2] /= mass_dm;
	vs_dm[0] /= mass_dm;
	vs_dm[1] /= mass_dm;
	vs_dm[2] /= mass_dm;
}
else
{
	s_dm[0] = Nodes[no].center[0];
	s_dm[1] = Nodes[no].center[1];
	s_dm[2] = Nodes[no].center[2];
	vs_dm[0] = 0;
	vs_dm[1] = 0;
	vs_dm[2] = 0;
}
#endif


Nodes[no].u.d.s[0] = s[0];
Nodes[no].u.d.s[1] = s[1];
Nodes[no].u.d.s[2] = s[2];
Extnodes[no].vs[0] = vs[0];
Extnodes[no].vs[1] = vs[1];
Extnodes[no].vs[2] = vs[2];
Nodes[no].u.d.mass = mass;
#ifdef RADTRANSFER
Nodes[no].stellar_s[0] = stellar_s[0];
Nodes[no].stellar_s[1] = stellar_s[1];
Nodes[no].stellar_s[2] = stellar_s[2];
Nodes[no].stellar_mass = stellar_mass;
#endif
#ifdef SCALARFIELD
Nodes[no].s_dm[0] = s_dm[0];
Nodes[no].s_dm[1] = s_dm[1];
Nodes[no].s_dm[2] = s_dm[2];
Nodes[no].mass_dm = mass_dm;
Extnodes[no].vs_dm[0] = vs_dm[0];
Extnodes[no].vs_dm[1] = vs_dm[1];
Extnodes[no].vs_dm[2] = vs_dm[2];
#endif

Extnodes[no].hmax = hmax;
Extnodes[no].vmax = vmax;
#ifndef KD_HMAX_ESTIMATE
Extnodes[no].divVmax = divVmax;
#endif

Extnodes[no].Flag = GlobFlag;

#ifdef FLTROUNDOFFREDUCTION
Extnodes[no].s_base[0] = Nodes[no].u.d.s[0];
Extnodes[no].s_base[1] = Nodes[no].u.d.s[1];
Extnodes[no].s_base[2] = Nodes[no].u.d.s[2];
#ifdef SCALARFIELD
Extnodes[no].s_dm_base[0] = Nodes[no].s_dm[0];
Extnodes[no].s_dm_base[1] = Nodes[no].s_dm[1];
Extnodes[no].s_dm_base[2] = Nodes[no].s_dm[2];
#endif
#endif

if(count_particles > 1)
	multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
else
	multiple_flag = 0;

Nodes[no].u.d.bitflags &= (~BITFLAG_MASK);	/* this clears the bits */

Nodes[no].u.d.bitflags |= multiple_flag;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
Nodes[no].u.d.bitflags |= diffsoftflag + (maxsofttype << BITFLAG_MAX_SOFTENING_TYPE);
#else
Nodes[no].maxsoft = maxsoft;
#endif
#endif

#ifdef ADAPTGRAVSOFT
Nodes[no].maxsoft = maxsoft;
#endif
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
	int no, i, m;

	/* mark all top-level nodes */

	for(i = 0; i < NTopleaves; i++)
	{
		no = DomainNodeIndex[i];

		while(no >= 0)
		{
			if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
				break;

			Nodes[no].u.d.bitflags |= (1 << BITFLAG_TOPLEVEL);

			no = Nodes[no].u.d.father;
		}

		/* mark also internal top level nodes */

		no = DomainNodeIndex[i];
		no = Nodes[no].u.d.father;

		while(no >= 0)
		{
			if(Nodes[no].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
				break;

			Nodes[no].u.d.bitflags |= (1 << BITFLAG_INTERNAL_TOPLEVEL);

			no = Nodes[no].u.d.father;
		}
	}

	/* mark top-level nodes that contain local particles */

	for(m = 0; m < MULTIPLEDOMAINS; m++)
		for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
				i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
		{
			no = DomainNodeIndex[i];

			if(DomainTask[i] != ThisTask)
				endrun(131231231);

			while(no >= 0)
			{
				if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))
					break;

				Nodes[no].u.d.bitflags |= (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS);

				no = Nodes[no].u.d.father;
			}
		}
}


/*! When a new additional star particle is created, we can put it into the
 *  tree at the position of the spawning gas particle. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed.
 */
void force_add_star_to_tree(int igas, int istar)
{
	int no;

	no = Nextnode[igas];
	Nextnode[igas] = istar;
	Nextnode[istar] = no;
	Father[istar] = Father[igas];
}






/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int force_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)
{
	struct NODE *nop = 0;
	int no, nexp, nodesinlist, ninteractions, ptype, task, listindex = 0;
	double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
	double pos_x, pos_y, pos_z, aold;
	MyLongDouble acc_x, acc_y, acc_z;

#ifdef SCALARFIELD
	double dx_dm = 0, dy_dm = 0, dz_dm = 0, mass_dm = 0;
#endif
#ifdef EVALPOTENTIAL
	double wp;
	MyLongDouble pot;

	pot = 0.0;
#endif

#ifdef FS_TURB_ESTIM
	int iturb;
	double Vel[4];
	MyDouble StrFnc[FS_BINS];
	long int StrFnc_count[FS_BINS];
#endif
#ifdef DISTORTIONTENSORPS
	int i1, i2;
	double fac2, h_tidal, h_inv_tidal, h3_inv_tidal, h5_inv_tidal;
	MyDouble tidal_tensorps[3][3];
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	double soft = 0;
#endif

#ifdef ADAPTGRAVSOFT
	double h_p, h_p_inv, h_p3_inv, u_p, fac_p, h_max;
	double zeta, omega, dWdr, dWdr_p, corr;
	double mass_target;
	int particle;
#ifdef AGS_OUTPUTCORR
	double correction = 0;
#endif
#endif
#ifdef DISTORTIONTENSORPS
	for(i1 = 0; i1 < 3; i1++)
		for(i2 = 0; i2 < 3; i2++)
			tidal_tensorps[i1][i2] = 0.0;
#endif

#ifdef FS_TURB_ESTIM
	for(iturb = 0; iturb < FS_BINS; iturb++)
	{
		StrFnc[iturb] = 0.;
		StrFnc_count[iturb] = 0;
	}
#endif

	acc_x = 0;
	acc_y = 0;
	acc_z = 0;
	ninteractions = 0;
	nodesinlist = 0;

	if(mode == 0)
	{
#ifdef GRAVITY_CENTROID
		if(P[target].Type == 0)
		{
			pos_x = SphP[target].Center[0];
			pos_y = SphP[target].Center[1];
			pos_z = SphP[target].Center[2];

		}
		else
		{
			pos_x = P[target].Pos[0];
			pos_y = P[target].Pos[1];
			pos_z = P[target].Pos[2];
		}
#else
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
#endif
#ifdef FS_TURB_ESTIM
		for(iturb = 0; iturb < 3; iturb++)
			Vel[iturb] = P[target].Vel[iturb];
		if(P[target].Type == 0)
			Vel[3] = PPP[target].Hsml;
		else
			Vel[3] = -1;
#endif

		ptype = P[target].Type;
		aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = PPP[target].Hsml;
#endif
#ifdef ADAPTGRAVSOFT
		zeta = P[target].AGS_zeta;
		omega = P[target].AGS_omega;
		h = P[target].AGS_Hsml;
		mass_target = P[target].Mass;
#endif
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#ifdef FS_TURB_ESTIM
		for(iturb = 0; iturb < 4; iturb++)
			Vel[iturb] = GravDataGet[target].Vel[iturb];
#endif
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		ptype = GravDataGet[target].Type;
#else
		ptype = P[0].Type;
#endif
		aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = GravDataGet[target].Soft;
#endif
#ifdef ADAPTGRAVSOFT
		zeta = GravDataGet[target].AGS_zeta;
		omega = GravDataGet[target].AGS_omega;
		h = GravDataGet[target].AGS_Hsml;
		mass_target = GravDataGet[target].Mass;
#endif
	}

#if !defined(UNEQUALSOFTENINGS) && !defined(ADAPTGRAVSOFT)
	h = All.ForceSoftening[ptype];
	h_inv = 1.0 / h;
	h3_inv = h_inv * h_inv * h_inv;
#endif

#ifdef DISTORTIONTENSORPS
	/* different tidal field softening */
	h_tidal = All.ForceSoftening[ptype];
	h_inv_tidal = 1.0 / h_tidal;
	h3_inv_tidal = h_inv_tidal * h_inv_tidal * h_inv_tidal;
	h5_inv_tidal = h_inv_tidal * h_inv_tidal * h_inv_tidal * h_inv_tidal * h_inv_tidal;
#endif


	if(mode == 0)
	{
		no = All.MaxPart;		/* root node */
	}
	else
	{
		nodesinlist++;
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < All.MaxPart)	/* single particle */
			{
				/* the index of the node is the index of the particle */
				/* observe the sign */

				if(P[no].Ti_current != All.Ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					drift_particle(no, All.Ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

#ifdef GRAVITY_CENTROID
				if(P[no].Type == 0)
				{
					dx = SphP[no].Center[0] - pos_x;
					dy = SphP[no].Center[1] - pos_y;
					dz = SphP[no].Center[2] - pos_z;

				}
				else
				{
					dx = P[no].Pos[0] - pos_x;
					dy = P[no].Pos[1] - pos_y;
					dz = P[no].Pos[2] - pos_z;
				}
#else
				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;
#endif
				mass = P[no].Mass;

#ifdef FS_TURB_ESTIM
				r = sqrt(dx * dx + dy * dy + dz * dz);
				if(Vel[3] > 0 && P[no].Type == 0 && r > 0)
				{
					iturb = (int) ceil(r / Vel[3] - 1.);
					if(iturb < FS_BINS)
					{
						StrFnc[iturb] +=
								fabs((Vel[0] * Vel[0] + Vel[1] * Vel[1] + Vel[2] * Vel[2]) -
										(P[no].Vel[0] * P[no].Vel[0] + P[no].Vel[1] * P[no].Vel[1] +
												P[no].Vel[2] * P[no].Vel[2]));
						StrFnc_count[iturb]++;
					}
				}
#endif

#ifdef SCALARFIELD
				if(ptype != 0)	/* we have a dark matter particle as target */
				{
					if(P[no].Type != 0)
					{
						dx_dm = dx;
						dy_dm = dy;
						dz_dm = dz;

						mass_dm = mass;
					}
					else
					{
						mass_dm = 0;
						dx_dm = dy_dm = dz_dm = 0;
					}
				}
#endif
			}
			else
			{
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							exportflag[task] = target;
							exportnodecount[task] = NODELISTLENGTH;
						}

						if(exportnodecount[task] == NODELISTLENGTH)
						{
							int exitFlag = 0;
							LOCK_NEXPORT;
#pragma omp critical(_nexport_)
							{
								if(Nexport >= All.BunchSize)
								{
									/* out of buffer space. Need to discard work for this particle and interrupt */
									BufferFullFlag = 1;
									exitFlag = 1;
								}
								else
								{
									nexp = Nexport;
									Nexport++;
								}
							}
							UNLOCK_NEXPORT;
							if(exitFlag)
								return -1;

							exportnodecount[task] = 0;
							exportindex[task] = nexp;
							DataIndexTable[nexp].Task = task;
							DataIndexTable[nexp].Index = target;
							DataIndexTable[nexp].IndexGet = nexp;
						}

						DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
								DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

						if(exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}

				nop = &Nodes[no];

				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				mass = nop->u.d.mass;

				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					if(mass)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}

				if(nop->Ti_current != All.Ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					force_drift_node(no, All.Ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;


#ifdef SCALARFIELD
				if(ptype != 0)	/* we have a dark matter particle as target */
				{
					dx_dm = nop->s_dm[0] - pos_x;
					dy_dm = nop->s_dm[1] - pos_y;
					dz_dm = nop->s_dm[2] - pos_z;

					mass_dm = nop->mass_dm;
				}
				else
				{
					mass_dm = 0;
					dx_dm = dy_dm = dz_dm = 0;
				}
#endif
			}

#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			dx = NEAREST(dx);
			dy = NEAREST(dy);
			dz = NEAREST(dz);
#endif
			r2 = dx * dx + dy * dy + dz * dz;


			if(no < All.MaxPart)
			{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];

				if(P[no].Type == 0)
				{
					if(h < PPP[no].Hsml)
						h = PPP[no].Hsml;
				}
				else
				{
					if(h < All.ForceSoftening[P[no].Type])
						h = All.ForceSoftening[P[no].Type];
				}
#else
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[P[no].Type])
					h = All.ForceSoftening[P[no].Type];
#endif
#endif
#ifdef ADAPTGRAVSOFT
				/* need to keep track of the softening (h_p) of the particle */
				h_p = P[no].AGS_Hsml;
				particle = no;
#endif
				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;
#ifdef THREAD_SAFE_COSTS
#pragma omp critical(_workcount_)
#endif
					P[no].GravCost[TakeLevel] += 1.0;
					UNLOCK_WORKCOUNT;
				}

				no = Nextnode[no];

#ifdef ADAPTGRAVSOFT
				/* otherwise r=0 and corr=nan. This check must be put after the next node has been selected! */
				if(target == particle)
					continue;
#endif
			}
			else			/* we have an  internal node. Need to check opening criterion */
			{

				if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
				{
					if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
				}
				else		/* check relative opening criterion */
				{
					if(mass * nop->len * nop->len > r2 * r2 * aold)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}

					/* check in addition whether we lie inside the cell */
					if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
					{
						if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
						{
							if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
							{
								no = nop->u.d.nextnode;
								continue;
							}
						}
					}
				}

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)])
				{
					h = All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)];
					if(r2 < h * h)
					{
						if(maskout_different_softening_flag(nop->u.d.bitflags))	/* signals that there are particles of different softening in the node */
						{
							no = nop->u.d.nextnode;
							continue;
						}
					}
				}
#else
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];

				if(h < nop->maxsoft)
				{
					h = nop->maxsoft;
					if(r2 < h * h)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}
#endif
#endif

#ifdef ADAPTGRAVSOFT
				h_p = nop->maxsoft;	/* softening of the node */
				h_max = h >= h_p ? h : h_p;

				if(r2 < h_max * h_max)
				{
					no = nop->u.d.nextnode;
					continue;	/* discard cases where the particle-node separation is less than the bigger
				   between the target softening and the node softening
					 */
				}

				particle = no;
#endif

				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;
#ifdef THREAD_SAFE_COSTS
#pragma omp critical(_workcount_)
#endif
					nop->GravCost += 1.0;
					UNLOCK_WORKCOUNT;
				}

				no = nop->u.d.sibling;	/* ok, node can be used */
			}

			r = sqrt(r2);

#ifdef ADAPTGRAVSOFT
			h_max = h >= h_p ? h : h_p;
			if(r >= h_max)
			{
				/* no need to worry about softening lengths: interaction is newtonian.
				 * All particle - node interactions fall here */
				fac = mass / (r2 * r);
				corr = 0;
			}
#else
			if(r >= h)
			{
				fac = mass / (r2 * r);
#ifdef EVALPOTENTIAL
				pot += FLT(-mass / r);
#endif
			}
#endif
			else
			{
#ifdef UNEQUALSOFTENINGS
				h_inv = 1.0 / h;
				h3_inv = h_inv * h_inv * h_inv;
#endif
#ifndef ADAPTGRAVSOFT
				u = r * h_inv;
				if(u < 0.5)
					fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
				else
					fac =
							mass * h3_inv * (21.333333333333 - 48.0 * u +
									38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
#ifdef EVALPOTENTIAL
				if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
				pot += FLT(mass * h_inv * wp);
#endif
#else
				/* interaction is smoothed: only particle-particle, no nodes down here! */
				dWdr = dWdr_p = corr = 0;
				h_inv = 1.0 / h;
				h3_inv = h_inv * h_inv * h_inv;
				u = r * h_inv;

				if(u > 1.)
					fac = mass / (r2 * r);

				else if(u < 0.5)
				{
					fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
					dWdr = 6. * KERNEL_COEFF_1 * h3_inv * h_inv * (-2. * u + 3. * u * u);
				}
				else
				{
					fac = mass * h3_inv * (21.333333333333 - 48.0 * u +
							38.4 * u * u - 10.666666666667 * u * u * u -
							0.066666666667 / (u * u * u));
					dWdr = -6. * KERNEL_COEFF_1 * h3_inv * h_inv * ((1. - u) * (1. - u));
				}


				h_p_inv = 1.0 / h_p;
				h_p3_inv = h_p_inv * h_p_inv * h_p_inv;
				u_p = r * h_p_inv;

				if(u_p > 1.)
					fac_p = mass / (r2 * r);

				else if(u_p < 0.5)
				{
					fac_p = mass * h_p3_inv * (10.666666666667 + u_p * u_p * (32.0 * u_p - 38.4));
					dWdr_p = 6. * KERNEL_COEFF_1 * h_p3_inv * h_p_inv * (-2. * u_p + 3. * u_p * u_p);
				}
				else
				{
					fac_p = mass * h_p3_inv * (21.333333333333 - 48.0 * u_p +
							38.4 * u_p * u_p - 10.666666666667 * u_p * u_p * u_p -
							0.066666666667 / (u_p * u_p * u_p));
					dWdr_p = -6. * KERNEL_COEFF_1 * h_p3_inv * h_p_inv * (1. - u_p) * (1. - u_p);
				}


				/* force = (fac + fac_p)/2 */
				fac += fac_p;
				fac /= 2.;

#ifndef AGS_NOCORRECTION
				/* correction term */
				corr =
						(zeta * omega * dWdr +
								mass / mass_target * P[particle].AGS_zeta * P[particle].AGS_omega * dWdr_p) / r;

				corr /= 2.;

#ifdef AGS_OUTPUTCORR
				correction += corr;
#endif
				/***************************************************************************************************/
				/* The original correction term is:                                                                */
				/*   corr = mass * (termI / termII * dWdr +  P[particle].termI / P[particle].termII * dWdr_p) / r; */
				/*   corr /= 2.;                                                                                   */
				/* But remember that termI and termII would have a different definition.                           */
				/***************************************************************************************************/
#endif
				if(particle >= All.MaxPart)
				{
					printf("*PROBLEM in ADAPTGRAVSOFT*: we are having a smoothed particle-node interaction!\n");
					printf("Particle, Node: %i, %i \n", target, particle);
				}
#endif
			}

#ifdef EVALPOTENTIAL
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
#endif

			acc_x += FLT(dx * fac);
			acc_y += FLT(dy * fac);
			acc_z += FLT(dz * fac);

#ifdef ADAPTGRAVSOFT
			acc_x += FLT(dx * corr);
			acc_y += FLT(dy * corr);
			acc_z += FLT(dz * corr);
#endif

			if(mass > 0)
				ninteractions++;

#ifdef SCALARFIELD
			if(ptype != 0)	/* we have a dark matter particle as target */
			{
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
				dx_dm = NEAREST(dx_dm);
				dy_dm = NEAREST(dy_dm);
				dz_dm = NEAREST(dz_dm);
#endif
				r2 = dx_dm * dx_dm + dy_dm * dy_dm + dz_dm * dz_dm;

				r = sqrt(r2);

				if(r >= h)
					fac = mass_dm / (r2 * r);
				else
				{
#ifdef UNEQUALSOFTENINGS
					h_inv = 1.0 / h;
					h3_inv = h_inv * h_inv * h_inv;
#endif
					u = r * h_inv;

					if(u < 0.5)
						fac = mass_dm * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
					else
						fac =
								mass_dm * h3_inv * (21.333333333333 - 48.0 * u +
										38.4 * u * u - 10.666666666667 * u * u * u -
										0.066666666667 / (u * u * u));
				}

				/* assemble force with strength, screening length, and target charge.  */

				fac *=
						All.ScalarBeta * (1 + r / All.ScalarScreeningLength) * exp(-r / All.ScalarScreeningLength);

				acc_x += FLT(dx_dm * fac);
				acc_y += FLT(dy_dm * fac);
				acc_z += FLT(dz_dm * fac);
			}
#endif


#ifdef DISTORTIONTENSORPS
			if(r >= h_tidal)
			{
				fac = mass / (r2 * r);
				fac2 = 3.0 * mass / (r2 * r2 * r);
			}
			else
			{
				u = r * h_inv_tidal;
				if(u < 0.5)
					fac = mass * h3_inv_tidal * (10.666666666667 + u * u * (32.0 * u - 38.4));
				else
					fac =
							mass * h3_inv_tidal * (21.333333333333 - 48.0 * u +
									38.4 * u * u - 10.666666666667 * u * u * u -
									0.066666666667 / (u * u * u));

				/*second derivates (see Gadget 1 paper and there g2 function) */
				if(u < 0.5)
					fac2 = mass * h5_inv_tidal * (76.8 - 96.0 * u);
				else
					fac2 = mass * h5_inv_tidal * (-0.2 / (u * u * u * u * u) + 48.0 / u - 76.8 + 32.0 * u);
			}


			/* tidal tensor */
			tidal_tensorps[0][0] += (-fac + dx * dx * fac2);
			tidal_tensorps[0][1] += (dx * dy * fac2);
			tidal_tensorps[0][2] += (dx * dz * fac2);
			tidal_tensorps[1][1] += (-fac + dy * dy * fac2);
			tidal_tensorps[1][2] += (dy * dz * fac2);
			tidal_tensorps[2][2] += (-fac + dz * dz * fac2);
			tidal_tensorps[1][0] = tidal_tensorps[0][1];
			tidal_tensorps[2][0] = tidal_tensorps[0][2];
			tidal_tensorps[2][1] = tidal_tensorps[1][2];
#endif



		}
		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
				{
					nodesinlist++;
					no = Nodes[no].u.d.nextnode;	/* open it */
				}
			}
		}
	}


	/* store result at the proper place */
	if(mode == 0)
	{
		P[target].g.dGravAccel[0] = acc_x;
		P[target].g.dGravAccel[1] = acc_y;
		P[target].g.dGravAccel[2] = acc_z;
#ifdef EVALPOTENTIAL
		P[target].p.dPotential = pot;
#endif
#ifdef FS_TURB_ESTIM
		for(iturb = 0; iturb < FS_BINS; iturb++)
		{
			P[target].StrFnc[iturb] = StrFnc[iturb];
			P[target].StrFnc_count[iturb] = StrFnc_count[iturb];
		}
#endif
#ifdef DISTORTIONTENSORPS
		for(i1 = 0; i1 < 3; i1++)
			for(i2 = 0; i2 < 3; i2++)
				P[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
		P[target].AGS_corr = correction;
#endif
	}
	else
	{
		GravDataResult[target].Acc[0] = acc_x;
		GravDataResult[target].Acc[1] = acc_y;
		GravDataResult[target].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
		GravDataResult[target].Potential = pot;
#endif
#ifdef FS_TURB_ESTIM
		for(iturb = 0; iturb < FS_BINS; iturb++)
		{
			GravDataResult[target].StrFnc[iturb] = StrFnc[iturb];
			GravDataResult[target].StrFnc_count[iturb] = StrFnc_count[iturb];
		}
#endif

#ifdef DISTORTIONTENSORPS
		for(i1 = 0; i1 < 3; i1++)
			for(i2 = 0; i2 < 3; i2++)
				GravDataResult[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
		GravDataResult[target].AGS_corr = correction;
#endif
		*exportflag = nodesinlist;
	}
	return ninteractions;
}

#ifdef PMGRID
/*! In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access panelty (which reduces cache performance) incurred by the
 *  table.
 */
int force_treeevaluate_shortrange(int target, int mode, int *exportflag, int *exportnodecount,
		int *exportindex)
{
	struct NODE *nop = 0;
	int no, nodesinlist, ptype, ninteractions, nexp, tabindex, task, listindex = 0;
	double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
	double dxx, dyy, dzz, pdxx, pdyy, pdzz;
	double pos_x, pos_y, pos_z, aold;
	double eff_dist;
	double rcut, asmth, asmthfac, rcut2, dist;
	MyLongDouble acc_x, acc_y, acc_z;
	// cache some global vars in local vars to help compiler with alias analysis
	int maxPart = All.MaxPart;
	int bunchSize = All.BunchSize;
	int maxNodes = MaxNodes;
	integertime ti_Current = All.Ti_Current;
	double errTol2 = All.ErrTolTheta * All.ErrTolTheta;

#ifdef PERIODIC
	1ST ////////////////////////////////////////////////////////////////////////////////////////////
	double xtmp;
#endif
#ifdef DISTORTIONTENSORPS
	2ND ////////////////////////////////////////////////////////////////////////////////////////////
	int i1, i2;
	double fac2, h5_inv;
	double fac_tidal;
	MyDouble tidal_tensorps[3][3];
#endif

#ifdef SCALARFIELD
	3RD ////////////////////////////////////////////////////////////////////////////////////////////
	double dx_dm = 0, dy_dm = 0, dz_dm = 0, mass_dm = 0;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	4TH ////////////////////////////////////////////////////////////////////////////////////////////
	double soft = 0;
#endif
#ifdef EVALPOTENTIAL
	5TH ////////////////////////////////////////////////////////////////////////////////////////////
	double wp, facpot;
	MyLongDouble pot;

	pot = 0;
#endif

#ifdef ADAPTGRAVSOFT
	6TH ////////////////////////////////////////////////////////////////////////////////////////////
	double h_p, h_p_inv, h_p3_inv, u_p, fac_p, h_max;
	double zeta, omega, dWdr, dWdr_p, corr;
	double mass_target;
	int particle;
#ifdef AGS_OUTPUTCORR
	7TH ////////////////////////////////////////////////////////////////////////////////////////////
	double correction = 0;
#endif
#endif

#ifdef DISTORTIONTENSORPS
	8TH ////////////////////////////////////////////////////////////////////////////////////////////
	for(i1 = 0; i1 < 3; i1++)
		for(i2 = 0; i2 < 3; i2++)
			tidal_tensorps[i1][i2] = 0.0;
#endif


	acc_x = 0;
	acc_y = 0;
	acc_z = 0;
	ninteractions = 0;
	nodesinlist = 0;

	rcut = All.Rcut[0];
	asmth = All.Asmth[0];

	if(mode != 0 && mode != 1)
	{
		printf("%d %d %d %d %d\n", target, mode, *exportflag, *exportnodecount, *exportindex);
		endrun(444);
	}

	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		ptype = P[target].Type;
		aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		9TH ////////////////////////////////////////////////////////////////////////////////////////////
		if(ptype == 0)
			soft = PPP[target].Hsml;
#endif
#ifdef PLACEHIGHRESREGION
		10TH ////////////////////////////////////////////////////////////////////////////////////////////
		if(pmforce_is_particle_high_res(ptype, P[target].Pos))
		{
			rcut = All.Rcut[1];
			asmth = All.Asmth[1];
		}
#endif
#ifdef ADAPTGRAVSOFT
		11TH ////////////////////////////////////////////////////////////////////////////////////////////
		zeta = P[target].AGS_zeta;
		omega = P[target].AGS_omega;
		h = P[target].AGS_Hsml;
		mass_target = P[target].Mass;
#endif
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		12TH ////////////////////////////////////////////////////////////////////////////////////////////
		ptype = GravDataGet[target].Type;
#else
		13TH ////////////////////////////////////////////////////////////////////////////////////////////
		ptype = P[0].Type;
#endif
		aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		14TH ////////////////////////////////////////////////////////////////////////////////////////////
		if(ptype == 0)
			soft = GravDataGet[target].Soft;
#endif
	}

	rcut2 = rcut * rcut;

	asmthfac = 0.5 / asmth * (NTAB / 3.0);

#if !defined(UNEQUALSOFTENINGS) && !defined(ADAPTGRAVSOFT)
	15TH ////////////////////////////////////////////////////////////////////////////////////////////
	h = All.ForceSoftening[ptype];
	h_inv = 1.0 / h;
	h3_inv = h_inv * h_inv * h_inv;
#ifdef DISTORTIONTENSORPS
	16TH ////////////////////////////////////////////////////////////////////////////////////////////
	h5_inv = h_inv * h_inv * h_inv * h_inv * h_inv;
#endif
#endif

	if(mode == 0)
	{
		no = maxPart;		/* root node */
	}
	else
	{
		nodesinlist++;
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < maxPart)
			{
				/* the index of the node is the index of the particle */
				if(P[no].Ti_current != ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					drift_particle(no, ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;

				dx = NEAREST(dx);
				dy = NEAREST(dy);
				dz = NEAREST(dz);

				r2 = dx * dx + dy * dy + dz * dz;

				mass = P[no].Mass;
#ifdef UNEQUALSOFTENINGS
				17TH ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
				18TH ////////////////////////////////////////////////////////////////////////////////////////////
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];

				if(P[no].Type == 0)
				{
					if(h < PPP[no].Hsml)
						h = PPP[no].Hsml;
				}
				else
				{
					if(h < All.ForceSoftening[P[no].Type])
						h = All.ForceSoftening[P[no].Type];
				}
#else
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[P[no].Type])
					h = All.ForceSoftening[P[no].Type];
#endif
#endif
				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;
#ifdef THREAD_SAFE_COSTS
					19TH ////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp critical(_workcount_)
#endif
					P[no].GravCost[TakeLevel] += 1.0;
					UNLOCK_WORKCOUNT;
				}
				no = Nextnode[no];
			}
			else			/* we have an  internal node */
			{
				if(no >= maxPart + maxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
						{
							exportflag[task] = target;
							exportnodecount[task] = NODELISTLENGTH;
						}

						if(exportnodecount[task] == NODELISTLENGTH)
						{
							int exitFlag = 0;
							LOCK_NEXPORT;
#pragma omp critical(_nexport_)
							{
								if(Nexport >= bunchSize)
								{
									/* out of buffer space. Need to discard work for this particle and interrupt */
									BufferFullFlag = 1;
									exitFlag = 1;
								}
								else
								{
									nexp = Nexport;
									Nexport++;
								}
							}
							UNLOCK_NEXPORT;
							if(exitFlag)
								return -1;

							exportnodecount[task] = 0;
							exportindex[task] = nexp;
							DataIndexTable[nexp].Task = task;
							DataIndexTable[nexp].Index = target;
							DataIndexTable[nexp].IndexGet = nexp;
						}

						DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
								DomainNodeIndex[no - (maxPart + maxNodes)];

						if(exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
					}
					no = Nextnode[no - maxNodes];
					continue;
				}

				nop = &Nodes[no];

				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					no = nop->u.d.nextnode;
					continue;
				}

				if(nop->Ti_current != ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					force_drift_node(no, ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				mass = nop->u.d.mass;

				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;

				dx = NEAREST(dx);
				dy = NEAREST(dy);
				dz = NEAREST(dz);
				r2 = dx * dx + dy * dy + dz * dz;
#ifdef DO_NOT_BRACH_IF
				20TH ////////////////////////////////////////////////////////////////////////////////////////////
				dxx = fabs(nop->center[0] - pos_x);
				dyy = fabs(nop->center[1] - pos_y);
				dzz = fabs(nop->center[2] - pos_z);
				eff_dist = rcut + 0.5 * nop->len;


				pdxx = NGB_PERIODIC_LONG_X(dxx);
				pdyy = NGB_PERIODIC_LONG_Y(dyy);
				pdzz = NGB_PERIODIC_LONG_Z(dzz);

				/* check whether we can stop walking along this branch */
				if((r2 > rcut2) & ((pdxx > eff_dist) | (pdyy > eff_dist) | (pdzz > eff_dist)))
				{
					no = nop->u.d.sibling;
					continue;
				}
#else
				21ST ////////////////////////////////////////////////////////////////////////////////////////////
				/* check whether we can stop walking along this branch */
				if(r2 > rcut2)
				{
					eff_dist = rcut + 0.5 * nop->len;

					dist = NEAREST(nop->center[0] - pos_x);
					if(dist < -eff_dist || dist > eff_dist)
					{
						no = nop->u.d.sibling;
						continue;
					}

					dist = NEAREST(nop->center[1] - pos_y);
					if(dist < -eff_dist || dist > eff_dist)
					{
						no = nop->u.d.sibling;
						continue;
					}

					dist = NEAREST(nop->center[2] - pos_z);
					if(dist < -eff_dist || dist > eff_dist)
					{
						no = nop->u.d.sibling;
						continue;
					}
				}
#endif

				if(errTol2)	/* check Barnes-Hut opening criterion */
				{
					if(nop->len * nop->len > r2 * errTol2)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
				}
				else		/* check relative opening criterion */
				{
#ifdef DO_NOT_BRACH_IF
					22ND ////////////////////////////////////////////////////////////////////////////////////////////
					if((mass * nop->len * nop->len > r2 * r2 * aold) |
							((dxx < 0.60 * nop->len) & (dyy < 0.60 * nop->len) & (dzz < 0.60 * nop->len)))
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
#else
					23RD ////////////////////////////////////////////////////////////////////////////////////////////
					if(mass * nop->len * nop->len > r2 * r2 * aold)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}

					/* check in addition whether we lie inside the cell */

					if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
					{
						if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
						{
							if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
							{
								no = nop->u.d.nextnode;
								continue;
							}
						}
					}
#endif
				}

#ifdef UNEQUALSOFTENINGS
				24TH ////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
				25TH ////////////////////////////////////////////////////////////////////////////////////////////
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)])
				{
					h = All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)];
					if(r2 < h * h)
					{
						if(maskout_different_softening_flag(nop->u.d.bitflags))	/* bit-5 signals that there are particles of different softening in the node */
						{
							no = nop->u.d.nextnode;

							continue;
						}
					}
				}
#else
				26TH ////////////////////////////////////////////////////////////////////////////////////////////
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];

				if(h < nop->maxsoft)
				{
					h = nop->maxsoft;
					if(r2 < h * h)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}
#endif
#endif
#ifdef ADAPTGRAVSOFT
				27TH ////////////////////////////////////////////////////////////////////////////////////////////
				h_p = nop->maxsoft;	/* softening of the node */
				h_max = h >= h_p ? h : h_p;

				if(r2 < h_max * h_max)
				{
					no = nop->u.d.nextnode;
					continue;	/* discard cases where the particle-node separation is less than the bigger
				   between the target softening and the node softening
					 */
				}

				particle = no;
#endif

				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;
#ifdef THREAD_SAFE_COSTS
					28TH ////////////////////////////////////////////////////////////////////////////////////////////
					29th  ////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp critical(_workcount_)
#endif
					nop->GravCost += 1.0;
					UNLOCK_WORKCOUNT;
				}

				no = nop->u.d.sibling;	/* ok, node can be used */

			}

			r = sqrt(r2);

#ifdef ADAPTGRAVSOFT
			30TH ////////////////////////////////////////////////////////////////////////////////////////////
			h_max = h >= h_p ? h : h_p;
			if(r >= h_max)
			{
				/* no need to worry about softening lengths: interaction is newtonian.
				 * All particle-node interactions fall here */
				fac = mass / (r2 * r);
				corr = 0;
			}
#else
			31ST ////////////////////////////////////////////////////////////////////////////////////////////
			if(r >= h)
			{
				fac = mass / (r2 * r);
#ifdef DISTORTIONTENSORPS
				32ND ////////////////////////////////////////////////////////////////////////////////////////////
				/* second derivative of potential needs this factor */
				fac2 = 3.0 * mass / (r2 * r2 * r);
#endif
#ifdef EVALPOTENTIAL
				33RD ////////////////////////////////////////////////////////////////////////////////////////////
				facpot = -mass / r;
#endif
			}
#endif
			else
			{
#ifndef ADAPTGRAVSOFT
				34TH  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef UNEQUALSOFTENINGS
				35TH  ////////////////////////////////////////////////////////////////////////////////////////////
				h_inv = 1.0 / h;
				h3_inv = h_inv * h_inv * h_inv;
#ifdef DISTORTIONTENSORPS
				36TH  ////////////////////////////////////////////////////////////////////////////////////////////
				h5_inv = h_inv * h_inv * h_inv * h_inv * h_inv;
#endif
#endif
				u = r * h_inv;
				if(u < 0.5)
					fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
				else
					fac =
							mass * h3_inv * (21.333333333333 - 48.0 * u +
									38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
#ifdef EVALPOTENTIAL
				37TH  ////////////////////////////////////////////////////////////////////////////////////////////
				if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

				facpot = mass * h_inv * wp;
#endif
#ifdef DISTORTIONTENSORPS
				38TH  ////////////////////////////////////////////////////////////////////////////////////////////
				/*second derivates needed -> calculate them from softend potential,
	         (see Gadget 1 paper and there g2 function). SIGN?! */
				if(u < 0.5)
					fac2 = mass * h5_inv * (76.8 - 96.0 * u);
				else
					fac2 = mass * h5_inv * (-0.2 / (u * u * u * u * u) + 48.0 / u - 76.8 + 32.0 * u);
#endif
#else
				39TH  ////////////////////////////////////////////////////////////////////////////////////////////
				/* interaction is smoothed: only particle-particle, no nodes down here! */
				dWdr = dWdr_p = corr = 0;
				h_inv = 1.0 / h;
				h3_inv = h_inv * h_inv * h_inv;
				u = r * h_inv;

				if(u > 1.)
					fac = mass / (r2 * r);

				else if(u < 0.5)
				{
					fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
					dWdr = 6. * KERNEL_COEFF_1 * h3_inv * h_inv * (-2. * u + 3. * u * u);
				}
				else
				{
					fac = mass * h3_inv * (21.333333333333 - 48.0 * u +
							38.4 * u * u - 10.666666666667 * u * u * u -
							0.066666666667 / (u * u * u));
					dWdr = -6. * KERNEL_COEFF_1 * h3_inv * h_inv * ((1. - u) * (1. - u));
				}


				h_p_inv = 1.0 / h_p;
				h_p3_inv = h_p_inv * h_p_inv * h_p_inv;
				u_p = r * h_p_inv;

				if(u_p > 1.)
					fac_p = mass / (r2 * r);

				else if(u_p < 0.5)
				{
					fac_p = mass * h_p3_inv * (10.666666666667 + u_p * u_p * (32.0 * u_p - 38.4));
					dWdr_p = 6. * KERNEL_COEFF_1 * h_p3_inv * h_p_inv * (-2. * u_p + 3. * u_p * u_p);
				}
				else
				{
					fac_p = mass * h_p3_inv * (21.333333333333 - 48.0 * u_p +
							38.4 * u_p * u_p - 10.666666666667 * u_p * u_p * u_p -
							0.066666666667 / (u_p * u_p * u_p));
					dWdr_p = -6. * KERNEL_COEFF_1 * h_p3_inv * h_p_inv * (1. - u_p) * (1. - u_p);
				}


				/* force = (fac + fac_p)/2 */
				fac += fac_p;
				fac /= 2.;

#ifndef AGS_NOCORRECTION
				40TH ////////////////////////////////////////////////////////////////////////////////////////////
				/* correction term */
				corr =
						(zeta * omega * dWdr +
								mass / mass_target * P[particle].AGS_zeta * P[particle].AGS_omega * dWdr_p) / r;

				corr /= 2.;


#ifdef AGS_OUTPUTCORR
				41ST  ////////////////////////////////////////////////////////////////////////////////////////////
				correction += corr;
#endif

				/***************************************************************************************************/
				/* The original correction term is:                                                                */
				/*   corr = mass * (termI / termII * dWdr +  P[particle].termI / P[particle].termII * dWdr_p) / r; */
				/*   corr /= 2.;                                                                                   */
				/* But remember that termI and termII would have a different definition.                           */
				/***************************************************************************************************/


#endif
				if(particle >= All.MaxPart)
				{
					printf("*PROBLEM in ADAPTGRAVSOFT*: we are having a smoothed particle-node interaction!\n");
					printf("Particle, Node: %i, %i \n", target, particle);
				}
#endif

			}

			tabindex = (int) (asmthfac * r);

			if(tabindex < NTAB)
			{
#ifdef DISTORTIONTENSORPS
				42ND  ////////////////////////////////////////////////////////////////////////////////////////////
				/* save original fac without shortrange_table facor (needed for tidal field calculation) */
				fac_tidal = fac;
#endif
				fac *= shortrange_table[tabindex];

				acc_x += FLT(dx * fac);
				acc_y += FLT(dy * fac);
				acc_z += FLT(dz * fac);

#ifdef ADAPTGRAVSOFT
				43RD  ////////////////////////////////////////////////////////////////////////////////////////////
				corr *= shortrange_table[tabindex];

				acc_x += FLT(dx * corr);
				acc_y += FLT(dy * corr);
				acc_z += FLT(dz * corr);
#endif

#ifdef DISTORTIONTENSORPS
				44TH  ////////////////////////////////////////////////////////////////////////////////////////////
				/*
	         tidal_tensorps[][] = Matrix of second derivatives of grav. potential, symmetric:
	         |Txx Txy Txz|   |tidal_tensorps[0][0] tidal_tensorps[0][1] tidal_tensorps[0][2]|
	         |Tyx Tyy Tyz| = |tidal_tensorps[1][0] tidal_tensorps[1][1] tidal_tensorps[1][2]| 
	         |Tzx Tzy Tzz|   |tidal_tensorps[2][0] tidal_tensorps[2][1] tidal_tensorps[2][2]|
				 */

				tidal_tensorps[0][0] += ((-fac_tidal + dx * dx * fac2) * shortrange_table[tabindex]) +
						dx * dx * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[0][1] += ((dx * dy * fac2) * shortrange_table[tabindex]) +
						dx * dy * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[0][2] += ((dx * dz * fac2) * shortrange_table[tabindex]) +
						dx * dz * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[1][1] += ((-fac_tidal + dy * dy * fac2) * shortrange_table[tabindex]) +
						dy * dy * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[1][2] += ((dy * dz * fac2) * shortrange_table[tabindex]) +
						dy * dz * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[2][2] += ((-fac_tidal + dz * dz * fac2) * shortrange_table[tabindex]) +
						dz * dz * fac2 / 3.0 * shortrange_table_tidal[tabindex];
				tidal_tensorps[1][0] = tidal_tensorps[0][1];
				tidal_tensorps[2][0] = tidal_tensorps[0][2];
				tidal_tensorps[2][1] = tidal_tensorps[1][2];
#endif
#ifdef EVALPOTENTIAL
				45TH  ////////////////////////////////////////////////////////////////////////////////////////////
				pot += FLT(facpot * shortrange_table_potential[tabindex]);
#endif
			}

			ninteractions++;


#ifdef SCALARFIELD
			46TH  ////////////////////////////////////////////////////////////////////////////////////////////
			if(ptype != 0)	/* we have a dark matter particle as target */
			{
#ifdef PERIODIC
				47TH  ////////////////////////////////////////////////////////////////////////////////////////////
				dx_dm = NEAREST(dx_dm);
				dy_dm = NEAREST(dy_dm);
				dz_dm = NEAREST(dz_dm);
#endif
				r2 = dx_dm * dx_dm + dy_dm * dy_dm + dz_dm * dz_dm;
				r = sqrt(r2);
				if(r >= h)
					fac = mass_dm / (r2 * r);
				else
				{
#ifdef UNEQUALSOFTENINGS
					48TH  ////////////////////////////////////////////////////////////////////////////////////////////
					h_inv = 1.0 / h;
					h3_inv = h_inv * h_inv * h_inv;
#endif
					u = r * h_inv;
					if(u < 0.5)
						fac = mass_dm * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
					else
						fac =
								mass_dm * h3_inv * (21.333333333333 - 48.0 * u +
										38.4 * u * u - 10.666666666667 * u * u * u -
										0.066666666667 / (u * u * u));
				}

				/* assemble force with strength, screening length, and target charge.  */

				fac *=
						All.ScalarBeta * (1 + r / All.ScalarScreeningLength) * exp(-r / All.ScalarScreeningLength);
				tabindex = (int) (asmthfac * r);
				if(tabindex < NTAB)
				{
					fac *= shortrange_table[tabindex];
					acc_x += FLT(dx_dm * fac);
					acc_y += FLT(dy_dm * fac);
					acc_z += FLT(dz_dm * fac);

				}
			}
#endif
		}

		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
				{
					nodesinlist++;
					no = Nodes[no].u.d.nextnode;	/* open it */
				}
			}
		}
	}


	/* store result at the proper place */
	if(mode == 0)
	{
		P[target].g.dGravAccel[0] = acc_x;
		P[target].g.dGravAccel[1] = acc_y;
		P[target].g.dGravAccel[2] = acc_z;
#ifdef EVALPOTENTIAL
		49TH  ////////////////////////////////////////////////////////////////////////////////////////////
		P[target].p.dPotential = pot;
#endif
#ifdef DISTORTIONTENSORPS
		50TH ////////////////////////////////////////////////////////////////////////////////////////////
		for(i1 = 0; i1 < 3; i1++)
			for(i2 = 0; i2 < 3; i2++)
				P[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
		51ST  ////////////////////////////////////////////////////////////////////////////////////////////
		P[target].AGS_corr = correction;
#endif
	}
	else
	{
		GravDataResult[target].Acc[0] = acc_x;
		GravDataResult[target].Acc[1] = acc_y;
		GravDataResult[target].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
		52ND  ////////////////////////////////////////////////////////////////////////////////////////////
		GravDataResult[target].Potential = pot;
#endif
#ifdef DISTORTIONTENSORPS
		53RD  ////////////////////////////////////////////////////////////////////////////////////////////
		for(i1 = 0; i1 < 3; i1++)
			for(i2 = 0; i2 < 3; i2++)
				GravDataResult[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];
#endif
#if defined(ADAPTGRAVSOFT) && defined(AGS_OUTPUTCORR)
		54TH  ////////////////////////////////////////////////////////////////////////////////////////////
		GravDataResult[target].AGS_corr = correction;
#endif
		*exportflag = nodesinlist;
	}

	return ninteractions;
}

#endif






#ifdef PERIODIC
/*! This function computes the Ewald correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly
 *  as it was done in Gadget-1.1. Instead, the tree is walked a second
 *  time. This is actually faster because the "Ewald-Treewalk" can use a
 *  different opening criterion than the normal tree walk. In particular,
 *  the Ewald correction is negligible for particles that are very close,
 *  but it is large for particles that are far away (this is quite
 *  different for the normal direct force). So we can here use a different
 *  opening criterion. Sufficient accuracy is usually obtained if the node
 *  length has dropped to a certain fraction ~< 0.25 of the
 *  BoxLength. However, we may only short-cut the interaction list of the
 *  normal full Ewald tree walk if we are sure that the whole node and all
 *  daughter nodes "lie on the same side" of the periodic boundary,
 *  i.e. that the real tree walk would not find a daughter node or particle
 *  that was mapped to a different nearest neighbour position when the tree
 *  walk would be further refined.
 */
int force_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount,
		int *exportindex)
{
	struct NODE *nop = 0;
	int no, cost, listindex = 0;
	double dx, dy, dz, mass, r2;
	int signx, signy, signz, nexp;
	int i, j, k, openflag, task;
	double u, v, w;
	double f1, f2, f3, f4, f5, f6, f7, f8;
	MyLongDouble acc_x, acc_y, acc_z;
	double boxsize, boxhalf;
	double pos_x, pos_y, pos_z, aold;

	boxsize = All.BoxSize;
	boxhalf = 0.5 * All.BoxSize;

	acc_x = 0;
	acc_y = 0;
	acc_z = 0;
	cost = 0;
	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		aold = All.ErrTolForceAcc * P[target].OldAcc;
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
		aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
	}

	if(mode == 0)
	{
		no = All.MaxPart;		/* root node */
	}
	else
	{
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < All.MaxPart)	/* single particle */
			{
				/* the index of the node is the index of the particle */
				/* observe the sign */
				if(P[no].Ti_current != All.Ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					drift_particle(no, All.Ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;
				mass = P[no].Mass;
			}
			else			/* we have an  internal node */
			{
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							exportflag[task] = target;
							exportnodecount[task] = NODELISTLENGTH;
						}

						if(exportnodecount[task] == NODELISTLENGTH)
						{
							int exitFlag = 0;
							LOCK_NEXPORT;
#pragma omp critical(_nexport_)
							{
								if(Nexport >= All.BunchSize)
								{
									/* out if buffer space. Need to discard work for this particle and interrupt */
									BufferFullFlag = 1;
									exitFlag = 1;
								}
								else
								{
									nexp = Nexport;
									Nexport++;
								}
							}
							UNLOCK_NEXPORT;
							if(exitFlag)
								return -1;

							exportnodecount[task] = 0;
							exportindex[task] = nexp;
							DataIndexTable[nexp].Task = task;
							DataIndexTable[nexp].Index = target;
							DataIndexTable[nexp].IndexGet = nexp;
						}

						DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
								DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

						if(exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}

				nop = &Nodes[no];

				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					no = nop->u.d.nextnode;
					continue;
				}

				if(nop->Ti_current != All.Ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					force_drift_node(no, All.Ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				mass = nop->u.d.mass;
				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;
			}

			dx = NEAREST(dx);
			dy = NEAREST(dy);
			dz = NEAREST(dz);

			if(no < All.MaxPart)
				no = Nextnode[no];
			else			/* we have an  internal node. Need to check opening criterion */
			{
				openflag = 0;
				r2 = dx * dx + dy * dy + dz * dz;
				if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
				{
					if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
					{
						openflag = 1;
					}
				}
				else		/* check relative opening criterion */
				{
					if(mass * nop->len * nop->len > r2 * r2 * aold)
					{
						openflag = 1;
					}
					else
					{
						if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
						{
							if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
							{
								if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
								{
									openflag = 1;
								}
							}
						}
					}
				}

				if(openflag)
				{
					/* now we check if we can avoid opening the cell */

					u = nop->center[0] - pos_x;
					if(u > boxhalf)
						u -= boxsize;
					if(u < -boxhalf)
						u += boxsize;
					if(fabs(u) > 0.5 * (boxsize - nop->len))
					{
						no = nop->u.d.nextnode;
						continue;
					}

					u = nop->center[1] - pos_y;
					if(u > boxhalf)
						u -= boxsize;
					if(u < -boxhalf)
						u += boxsize;
					if(fabs(u) > 0.5 * (boxsize - nop->len))
					{
						no = nop->u.d.nextnode;
						continue;
					}

					u = nop->center[2] - pos_z;
					if(u > boxhalf)
						u -= boxsize;
					if(u < -boxhalf)
						u += boxsize;
					if(fabs(u) > 0.5 * (boxsize - nop->len))
					{
						no = nop->u.d.nextnode;
						continue;
					}

					/* if the cell is too large, we need to refine
					 * it further
					 */
					if(nop->len > 0.20 * boxsize)
					{
						/* cell is too large */
						no = nop->u.d.nextnode;
						continue;
					}
				}

				no = nop->u.d.sibling;	/* ok, node can be used */
			}

			/* compute the Ewald correction force */

			if(dx < 0)
			{
				dx = -dx;
				signx = +1;
			}
			else
				signx = -1;
			if(dy < 0)
			{
				dy = -dy;
				signy = +1;
			}
			else
				signy = -1;
			if(dz < 0)
			{
				dz = -dz;
				signz = +1;
			}
			else
				signz = -1;
			u = dx * fac_intp;
			i = (int) u;
			if(i >= EN)
				i = EN - 1;
			u -= i;
			v = dy * fac_intp;
			j = (int) v;
			if(j >= EN)
				j = EN - 1;
			v -= j;
			w = dz * fac_intp;
			k = (int) w;
			if(k >= EN)
				k = EN - 1;
			w -= k;
			/* compute factors for trilinear interpolation */
			f1 = (1 - u) * (1 - v) * (1 - w);
			f2 = (1 - u) * (1 - v) * (w);
			f3 = (1 - u) * (v) * (1 - w);
			f4 = (1 - u) * (v) * (w);
			f5 = (u) * (1 - v) * (1 - w);
			f6 = (u) * (1 - v) * (w);
			f7 = (u) * (v) * (1 - w);
			f8 = (u) * (v) * (w);
			acc_x += FLT(mass * signx * (fcorrx[i][j][k] * f1 +
					fcorrx[i][j][k + 1] * f2 +
					fcorrx[i][j + 1][k] * f3 +
					fcorrx[i][j + 1][k + 1] * f4 +
					fcorrx[i + 1][j][k] * f5 +
					fcorrx[i + 1][j][k + 1] * f6 +
					fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8));
			acc_y +=
					FLT(mass * signy *
							(fcorry[i][j][k] * f1 + fcorry[i][j][k + 1] * f2 +
									fcorry[i][j + 1][k] * f3 + fcorry[i][j + 1][k + 1] * f4 + fcorry[i +
									                                                                 1]
									                                                                 [j][k] * f5 + fcorry[i + 1][j][k + 1] * f6 + fcorry[i + 1][j +
									                                                                                                                            1][k] *
									                                                                                                                            f7 + fcorry[i + 1][j + 1][k + 1] * f8));
			acc_z +=
					FLT(mass * signz *
							(fcorrz[i][j][k] * f1 + fcorrz[i][j][k + 1] * f2 +
									fcorrz[i][j + 1][k] * f3 + fcorrz[i][j + 1][k + 1] * f4 + fcorrz[i +
									                                                                 1]
									                                                                 [j][k] * f5 + fcorrz[i + 1][j][k + 1] * f6 + fcorrz[i + 1][j +
									                                                                                                                            1][k] *
									                                                                                                                            f7 + fcorrz[i + 1][j + 1][k + 1] * f8));
			cost++;
		}

		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
					no = Nodes[no].u.d.nextnode;	/* open it */
			}
		}
	}

	/* add the result at the proper place */

	if(mode == 0)
	{
		P[target].g.dGravAccel[0] += acc_x;
		P[target].g.dGravAccel[1] += acc_y;
		P[target].g.dGravAccel[2] += acc_z;
	}
	else
	{
		GravDataResult[target].Acc[0] = acc_x;
		GravDataResult[target].Acc[1] = acc_y;
		GravDataResult[target].Acc[2] = acc_z;
	}

	return cost;
}

#endif






/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
int force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
	struct NODE *nop = 0;
	MyLongDouble pot;
	int no, ptype, task, nexport_save, listindex = 0;
	double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
	double pos_x, pos_y, pos_z, aold;

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	double soft = 0;
#endif

#ifdef ADAPTGRAVSOFT
	double h_p, h_p_inv, u_p, wp_p, h_max;
	int particle;
#endif

	nexport_save = *nexport;
	pot = 0;
	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		ptype = P[target].Type;
		aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = PPP[target].Hsml;
#endif
#ifdef ADAPTGRAVSOFT
		h = P[target].AGS_Hsml;
#endif
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		ptype = GravDataGet[target].Type;
#else
		ptype = P[0].Type;
#endif
		aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = GravDataGet[target].Soft;
#endif
#ifdef ADAPTGRAVSOFT
		h = GravDataGet[target].AGS_Hsml;
#endif
	}

#if !defined(UNEQUALSOFTENINGS) && !defined(ADAPTGRAVSOFT)
	h = All.ForceSoftening[ptype];
	h_inv = 1.0 / h;
#endif
	if(mode == 0)
	{
		no = All.MaxPart;		/* root node */
	}
	else
	{
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < All.MaxPart)	/* single particle */
			{
				/* the index of the node is the index of the particle */
				/* observe the sign */

				if(P[no].Ti_current != All.Ti_Current)
					drift_particle(no, All.Ti_Current);
				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;
				mass = P[no].Mass;
			}
			else
			{
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							Exportflag[task] = target;
							Exportnodecount[task] = NODELISTLENGTH;
						}

						if(Exportnodecount[task] == NODELISTLENGTH)
						{
							if(*nexport >= All.BunchSize)
							{
								*nexport = nexport_save;
								if(nexport_save == 0)
									endrun(13001);	/* in this case, the buffer is too small to process even a single particle */
								for(task = 0; task < NTask; task++)
									nsend_local[task] = 0;
								for(no = 0; no < nexport_save; no++)
									nsend_local[DataIndexTable[no].Task]++;
								return -1;
							}
							Exportnodecount[task] = 0;
							Exportindex[task] = *nexport;
							DataIndexTable[*nexport].Task = task;
							DataIndexTable[*nexport].Index = target;
							DataIndexTable[*nexport].IndexGet = *nexport;
							*nexport = *nexport + 1;
							nsend_local[task]++;
						}

						DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
								DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
						if(Exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}

				nop = &Nodes[no];
				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					no = nop->u.d.nextnode;
					continue;
				}

				if(nop->Ti_current != All.Ti_Current)
					force_drift_node(no, All.Ti_Current);
				mass = nop->u.d.mass;
				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;
			}

#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			dx = NEAREST(dx);
			dy = NEAREST(dy);
			dz = NEAREST(dz);
#endif
			r2 = dx * dx + dy * dy + dz * dz;
			if(no < All.MaxPart)
			{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];
				if(P[no].Type == 0)
				{
					if(h < PPP[no].Hsml)
						h = PPP[no].Hsml;
				}
				else
				{
					if(h < All.ForceSoftening[P[no].Type])
						h = All.ForceSoftening[P[no].Type];
				}
#else
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[P[no].Type])
					h = All.ForceSoftening[P[no].Type];
#endif
#endif
#ifdef ADAPTGRAVSOFT
				/* need to keep track of the softening (h_p) of the particle */
				h_p = P[no].AGS_Hsml;
				particle = no;
#endif
				no = Nextnode[no];
			}
			else			/* we have an internal node. Need to check opening criterion */
			{
				if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
				{
					if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
				}
				else		/* check relative opening criterion */
				{
					if(mass * nop->len * nop->len > r2 * r2 * aold)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}

					if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
					{
						if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
						{
							if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
							{
								no = nop->u.d.nextnode;
								continue;
							}
						}
					}
				}
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)])
				{
					h = All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)];
					if(r2 < h * h)
					{
						if(maskout_different_softening_flag(nop->u.d.bitflags))	/* bit-5 signals that there are particles of different softening in the node */
						{
							no = nop->u.d.nextnode;
							continue;
						}
					}
				}
#else
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];
				if(h < nop->maxsoft)
				{
					h = nop->maxsoft;
					if(r2 < h * h)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}
#endif
#endif
#ifdef ADAPTGRAVSOFT

				h_p = nop->maxsoft;	/* softening of the node */
				h_max = h >= h_p ? h : h_p;

				if(r2 < h_max * h_max)
				{
					no = nop->u.d.nextnode;
					continue;	/* discard cases where the particle-node separation is less than the bigger
				   between the target softening and the node softening
					 */
				}

				particle = no;
#endif
				no = nop->u.d.sibling;	/* node can be used */
			}

			r = sqrt(r2);

#ifdef ADAPTGRAVSOFT
			h_max = h >= h_p ? h : h_p;
			if(r >= h_max)
				pot += FLT(-mass / r);
#else
			if(r >= h)
				pot += FLT(-mass / r);
#endif
			else
			{
#ifdef UNEQUALSOFTENINGS
				h_inv = 1.0 / h;
#endif
#ifdef ADAPTGRAVSOFT

				h_inv = 1.0 / h;
				u = r * h_inv;

				if(u > 1.)
					wp = -1. / u;
				else if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

				h_p_inv = 1. / h_p;
				u_p = r * h_p_inv;

				if(u_p > 1.)
					wp_p = -1. / u_p;
				else if(u_p < 0.5)
					wp_p = -2.8 + u_p * u_p * (5.333333333333 + u_p * u_p * (6.4 * u_p - 9.6));
				else
					wp_p =
							-3.2 + 0.066666666667 / u_p + u_p * u_p * (10.666666666667 +
									u_p * (-16.0 +
											u_p * (9.6 - 2.133333333333 * u_p)));


				pot += FLT((mass * h_inv * wp + mass * h_p_inv * wp_p) / 2.);


#else
				u = r * h_inv;
				if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
				pot += FLT(mass * h_inv * wp);
#endif
			}
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
		}
		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
					no = Nodes[no].u.d.nextnode;	/* open it */
			}
		}
	}

	/* store result at the proper place */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
	if(mode == 0)
		P[target].p.dPotential = pot;
	else
		PotDataResult[target].Potential = pot;
#endif
	return 0;
}



#ifdef SUBFIND
int subfind_force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
	struct NODE *nop = 0;
	MyLongDouble pot;
	int no, ptype, task, nexport_save, listindex = 0;
	double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
	double pos_x, pos_y, pos_z;
#ifdef ADAPTGRAVSOFT
	double h_p, h_p_inv, u_p, wp_p, h_max;
	int particle;
#endif

	nexport_save = *nexport;
	pot = 0;
	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		ptype = P[target].Type;
#ifdef ADAPTGRAVSOFT
		h = P[target].AGS_Hsml;
#endif
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		ptype = GravDataGet[target].Type;
#else
		ptype = P[0].Type;
#endif
#ifdef ADAPTGRAVSOFT
		h = GravDataGet[target].AGS_Hsml;
#endif
	}

#ifndef ADAPTGRAVSOFT
	h = All.ForceSoftening[ptype];
	h_inv = 1.0 / h;
#endif

	if(mode == 0)
	{
		no = All.MaxPart;		/* root node */
	}
	else
	{
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < All.MaxPart)	/* single particle */
			{
				/* the index of the node is the index of the particle */
				/* observe the sign */

				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;
				mass = P[no].Mass;
			}
			else
			{
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							Exportflag[task] = target;
							Exportnodecount[task] = NODELISTLENGTH;
						}

						if(Exportnodecount[task] == NODELISTLENGTH)
						{
							if(*nexport >= All.BunchSize)
							{
								*nexport = nexport_save;
								if(nexport_save == 0)
									endrun(13001);	/* in this case, the buffer is too small to process even a single particle */
								for(task = 0; task < NTask; task++)
									nsend_local[task] = 0;
								for(no = 0; no < nexport_save; no++)
									nsend_local[DataIndexTable[no].Task]++;
								return -1;
							}
							Exportnodecount[task] = 0;
							Exportindex[task] = *nexport;
							DataIndexTable[*nexport].Task = task;
							DataIndexTable[*nexport].Index = target;
							DataIndexTable[*nexport].IndexGet = *nexport;
							*nexport = *nexport + 1;
							nsend_local[task]++;
						}

						DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
								DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
						if(Exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}

				nop = &Nodes[no];
				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				mass = nop->u.d.mass;
				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					if(mass)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}

				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;
			}

#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			dx = NEAREST(dx);
			dy = NEAREST(dy);
			dz = NEAREST(dz);
#endif
			r2 = dx * dx + dy * dy + dz * dz;
			if(no < All.MaxPart)
			{
#ifdef ADAPTGRAVSOFT
				/* need to keep track of the softening (h_p) of the particle */
				h_p = P[no].AGS_Hsml;
				particle = no;
#endif
				no = Nextnode[no];
			}
			else			/* we have an internal node. Need to check opening criterion */
			{
				/* check Barnes-Hut opening criterion */

				if(nop->len * nop->len > r2 * All.ErrTolThetaSubfind * All.ErrTolThetaSubfind)
				{
					/* open cell */
					if(mass)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}
#ifdef ADAPTGRAVSOFT
				h_p = nop->maxsoft;	/* softening of the node */
				h_max = h >= h_p ? h : h_p;

				if(r2 < h_max * h_max)
				{
					no = nop->u.d.nextnode;
					continue;	/* discard cases where the particle-node separation is less than the bigger
				   between the target softening and the node softening
					 */
				}

				particle = no;
#endif
				no = nop->u.d.sibling;	/* node can be used */
			}

			r = sqrt(r2);
#ifdef ADAPTGRAVSOFT
			h_max = h >= h_p ? h : h_p;
			if(r >= h_max)
				pot += FLT(-mass / r);
#else
			if(r >= h)
				pot += FLT(-mass / r);
#endif
			else
			{
#ifdef ADAPTGRAVSOFT
				h_inv = 1.0 / h;
				u = r * h_inv;

				if(u > 1.)
					wp = -1. / u;
				else if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

				h_p_inv = 1. / h_p;
				u_p = r * h_p_inv;

				if(u_p > 1.)
					wp_p = -1. / u_p;
				else if(u_p < 0.5)
					wp_p = -2.8 + u_p * u_p * (5.333333333333 + u_p * u_p * (6.4 * u_p - 9.6));
				else
					wp_p =
							-3.2 + 0.066666666667 / u_p + u_p * u_p * (10.666666666667 +
									u_p * (-16.0 +
											u_p * (9.6 - 2.133333333333 * u_p)));

				pot += FLT((mass * h_inv * wp + mass * h_p_inv * wp_p) / 2.);
#else
				u = r * h_inv;
				if(u < 0.5)
					wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
				else
					wp =
							-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
									u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
				pot += FLT(mass * h_inv * wp);
#endif
			}
		}
		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
					no = Nodes[no].u.d.nextnode;	/* open it */
			}
		}
	}

	/* store result at the proper place */

	if(mode == 0)
		P[target].u.DM_Potential = pot;
	else
		PotDataResult[target].Potential = pot;
	return 0;
}
#endif



#ifdef PMGRID
/*! This function computes the short-range potential when the TreePM
 *  algorithm is used. This potential is the Newtonian potential, modified
 *  by a complementary error function.
 */
int force_treeevaluate_potential_shortrange(int target, int mode, int *nexport, int *nsend_local)
{
	struct NODE *nop = 0;
	MyLongDouble pot;
	int no, ptype, tabindex, task, nexport_save, listindex = 0;
	double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
	double pos_x, pos_y, pos_z, aold;
	double eff_dist, fac, rcut, asmth, asmthfac;
	double dxx, dyy, dzz;

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	double soft = 0;
#endif
#ifdef ADAPTGRAVSOFT
	double h_p, h_p_inv, u_p, wp_p, h_max;
	int particle;
#endif

	nexport_save = *nexport;
	pot = 0;

	rcut = All.Rcut[0];
	asmth = All.Asmth[0];

	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		ptype = P[target].Type;
		aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = PPP[target].Hsml;
#endif
#ifdef PLACEHIGHRESREGION
		if(pmforce_is_particle_high_res(ptype, P[target].Pos))
		{
			rcut = All.Rcut[1];
			asmth = All.Asmth[1];
		}
#endif
#ifdef ADAPTGRAVSOFT
		h = P[target].AGS_Hsml;
#endif
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		ptype = GravDataGet[target].Type;
#else
		ptype = P[0].Type;
#endif
		aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		if(ptype == 0)
			soft = GravDataGet[target].Soft;
#endif
#ifdef PLACEHIGHRESREGION
		if(pmforce_is_particle_high_res(ptype, GravDataGet[target].Pos))
		{
			rcut = All.Rcut[1];
			asmth = All.Asmth[1];
		}
#endif
#ifdef ADAPTGRAVSOFT
		h = GravDataGet[target].AGS_Hsml;
#endif
	}


	asmthfac = 0.5 / asmth * (NTAB / 3.0);
#if !defined(UNEQUALSOFTENINGS) && !defined(ADAPTGRAVSOFT)
	h = All.ForceSoftening[ptype];
	h_inv = 1.0 / h;
#endif
	if(mode == 0)
	{
		no = All.MaxPart;		/* root node */
	}
	else
	{
		no = GravDataGet[target].NodeList[0];
		no = Nodes[no].u.d.nextnode;	/* open it */
	}

	while(no >= 0)
	{
		while(no >= 0)
		{
			if(no < All.MaxPart)	/* single particle */
			{
				/* the index of the node is the index of the particle */
				/* observe the sign  */
#ifndef SUBFIND_RESHUFFLE_AND_POTENTIAL
				if(P[no].Ti_current != All.Ti_Current)
					drift_particle(no, All.Ti_Current);
#endif
				dx = P[no].Pos[0] - pos_x;
				dy = P[no].Pos[1] - pos_y;
				dz = P[no].Pos[2] - pos_z;
				mass = P[no].Mass;
			}
			else
			{
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							Exportflag[task] = target;
							Exportnodecount[task] = NODELISTLENGTH;
						}

						if(Exportnodecount[task] == NODELISTLENGTH)
						{
							if(*nexport >= All.BunchSize)
							{
								*nexport = nexport_save;
								if(nexport_save == 0)
									endrun(13002);	/* in this case, the buffer is too small to process even a single particle */
								for(task = 0; task < NTask; task++)
									nsend_local[task] = 0;
								for(no = 0; no < nexport_save; no++)
									nsend_local[DataIndexTable[no].Task]++;
								return -1;
							}
							Exportnodecount[task] = 0;
							Exportindex[task] = *nexport;
							DataIndexTable[*nexport].Task = task;
							DataIndexTable[*nexport].Index = target;
							DataIndexTable[*nexport].IndexGet = *nexport;
							*nexport = *nexport + 1;
							nsend_local[task]++;
						}

						DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
								DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
						if(Exportnodecount[task] < NODELISTLENGTH)
							DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}

				nop = &Nodes[no];
				if(mode == 1)
				{
					if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						no = -1;
						continue;
					}
				}

				if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					no = nop->u.d.nextnode;
					continue;
				}
#ifndef SUBFIND_RESHUFFLE_AND_POTENTIAL
				if(nop->Ti_current != All.Ti_Current)
					force_drift_node(no, All.Ti_Current);
#endif
				mass = nop->u.d.mass;
				dx = nop->u.d.s[0] - pos_x;
				dy = nop->u.d.s[1] - pos_y;
				dz = nop->u.d.s[2] - pos_z;
			}

#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
			dx = NEAREST(dx);
			dy = NEAREST(dy);
			dz = NEAREST(dz);
#endif
			r2 = dx * dx + dy * dy + dz * dz;
			if(no < All.MaxPart)
			{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];
				if(P[no].Type == 0)
				{
					if(h < PPP[no].Hsml)
						h = PPP[no].Hsml;
				}
				else
				{
					if(h < All.ForceSoftening[P[no].Type])
						h = All.ForceSoftening[P[no].Type];
				}
#else
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[P[no].Type])
					h = All.ForceSoftening[P[no].Type];
#endif
#endif
#ifdef ADAPTGRAVSOFT
				/* need to keep track of the softening (h_p) of the particle */
				h_p = P[no].AGS_Hsml;
				particle = no;
#endif
				no = Nextnode[no];
			}
			else			/* we have an  internal node. Need to check opening criterion */
			{
				/* check whether we can stop walking along this branch */
				if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
				{
					if(mode == 0)
					{
						if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
						{
							Exportflag[task] = target;
							DataIndexTable[*nexport].Index = target;
							DataIndexTable[*nexport].Task = task;	/* Destination task */
							*nexport = *nexport + 1;
							nsend_local[task]++;
						}
					}
					no = Nextnode[no - MaxNodes];
					continue;
				}


				eff_dist = rcut + 0.5 * nop->len;
				dxx = nop->center[0] - pos_x;	/* observe the sign ! */
				dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
				dzz = nop->center[2] - pos_z;
#ifdef PERIODIC
				dxx = NEAREST(dxx);
				dyy = NEAREST(dyy);
				dzz = NEAREST(dzz);
#endif
#ifdef DO_NOT_BRACH_IF
				if((fabs(dxx) > eff_dist) | (fabs(dyy) > eff_dist) | (fabs(dzz) > eff_dist))
				{
					no = nop->u.d.sibling;
					continue;
				}
#else
				if(dxx < -eff_dist || dxx > eff_dist)
				{
					no = nop->u.d.sibling;
					continue;
				}

				if(dyy < -eff_dist || dyy > eff_dist)
				{
					no = nop->u.d.sibling;
					continue;
				}

				if(dzz < -eff_dist || dzz > eff_dist)
				{
					no = nop->u.d.sibling;
					continue;
				}
#endif

				if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
				{
					if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
				}
				else		/* check relative opening criterion */
				{
#ifdef DO_NOT_BRACH_IF
					if((mass * nop->len * nop->len > r2 * r2 * aold) |
							((fabs(dxx) < 0.60 * nop->len) & (fabs(dyy) < 0.60 * nop->len) & (fabs(dzz) <
									0.60 * nop->len)))
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}
#else
					if(mass * nop->len * nop->len > r2 * r2 * aold)
					{
						/* open cell */
						no = nop->u.d.nextnode;
						continue;
					}

					if(fabs(dxx) < 0.60 * nop->len)
					{
						if(fabs(dyy) < 0.60 * nop->len)
						{
							if(fabs(dzz) < 0.60 * nop->len)
							{
								no = nop->u.d.nextnode;
								continue;
							}
						}
					}
#endif
				}

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
				h = All.ForceSoftening[ptype];
				if(h < All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)])
				{
					h = All.ForceSoftening[extract_max_softening_type(nop->u.d.bitflags)];
					if(r2 < h * h)
					{
						/* bit-5 signals that there are particles of
						 * different softening in the node
						 */
						if(maskout_different_softening_flag(nop->u.d.bitflags))
						{
							no = nop->u.d.nextnode;
							continue;
						}
					}
				}
#else
				if(ptype == 0)
					h = soft;
				else
					h = All.ForceSoftening[ptype];
				if(h < nop->maxsoft)
				{
					h = nop->maxsoft;
					if(r2 < h * h)
					{
						no = nop->u.d.nextnode;
						continue;
					}
				}
#endif
#endif
#ifdef ADAPTGRAVSOFT

				h_p = nop->maxsoft;	/* softening of the node */
				h_max = h >= h_p ? h : h_p;

				if(r2 < h_max * h_max)
				{
					no = nop->u.d.nextnode;
					continue;	/* discard cases where the particle-node separation is less than the bigger
				   between the target softening and the node softening
					 */
				}

				particle = no;
#endif
				no = nop->u.d.sibling;	/* node can be used */
			}

			r = sqrt(r2);
			tabindex = (int) (r * asmthfac);
			if(tabindex < NTAB)
			{
				fac = shortrange_table_potential[tabindex];
#ifdef ADAPTGRAVSOFT
				h_max = h >= h_p ? h : h_p;
				if(r >= h_max)
					pot += FLT(-fac * mass / r);
#else
				if(r >= h)
					pot += FLT(-fac * mass / r);
#endif
				else
				{
#ifdef UNEQUALSOFTENINGS
					h_inv = 1.0 / h;
#endif
#ifdef ADAPTGRAVSOFT
					h_inv = 1.0 / h;
					u = r * h_inv;

					if(u > 1.)
						wp = -1. / u;
					else if(u < 0.5)
						wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
					else
						wp =
								-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
										u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

					h_p_inv = 1.0 / h_p;
					u_p = r * h_p_inv;

					if(u_p > 1.)
						wp_p = -1. / u_p;
					else if(u_p < 0.5)
						wp_p = -2.8 + u_p * u_p * (5.333333333333 + u_p * u_p * (6.4 * u_p - 9.6));
					else
						wp_p =
								-3.2 + 0.066666666667 / u_p + u_p * u_p * (10.666666666667 +
										u_p * (-16.0 +
												u_p * (9.6 - 2.133333333333 * u_p)));


					pot += FLT(fac * (mass * h_inv * wp + mass * h_p_inv * wp_p) / 2.);

#else
					u = r * h_inv;
					if(u < 0.5)
						wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
					else
						wp =
								-3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
										u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
					pot += FLT(fac * mass * h_inv * wp);
#endif
				}
			}
		}
		if(mode == 1)
		{
			listindex++;
			if(listindex < NODELISTLENGTH)
			{
				no = GravDataGet[target].NodeList[listindex];
				if(no >= 0)
					no = Nodes[no].u.d.nextnode;	/* open it */
			}
		}
	}

	/* store result at the proper place */
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
	if(mode == 0)
		P[target].p.dPotential = pot;
	else
		PotDataResult[target].Potential = pot;
#endif
	return 0;
}

#endif



/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
	int i;
	size_t bytes;
	double allbytes = 0, allbytes_topleaves = 0;
	double u;

	tree_allocated_flag = 1;
	DomainNodeIndex = (int *) mymalloc("DomainNodeIndex", bytes = NTopleaves * sizeof(int));
	allbytes_topleaves += bytes;
	MaxNodes = maxnodes;
	if(!(Nodes_base = (struct NODE *) mymalloc("Nodes_base", bytes = (MaxNodes + 1) * sizeof(struct NODE))))
	{
		printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
		endrun(3);
	}
	allbytes += bytes;
	if(!
			(Extnodes_base =
					(struct extNODE *) mymalloc("Extnodes_base", bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
	{
		printf("failed to allocate memory for %d tree-extnodes (%g MB).\n",
				MaxNodes, bytes / (1024.0 * 1024.0));
		endrun(3);
	}
	allbytes += bytes;
	Nodes = Nodes_base - All.MaxPart;
	Extnodes = Extnodes_base - All.MaxPart;
	if(!(Nextnode = (int *) mymalloc("Nextnode", bytes = (maxpart + NTopnodes) * sizeof(int))))
	{
		printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n",
				maxpart + NTopnodes, bytes / (1024.0 * 1024.0));
		exit(0);
	}
	allbytes += bytes;
	if(!(Father = (int *) mymalloc("Father", bytes = (maxpart) * sizeof(int))))
	{
		printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
		exit(0);
	}
	allbytes += bytes;
	if(first_flag == 0)
	{
		first_flag = 1;
		if(ThisTask == 0)
			printf
			("\nAllocated %g MByte for BH-tree, and %g Mbyte for top-leaves.  (presently allocted %g MB)\n\n",
					allbytes / (1024.0 * 1024.0), allbytes_topleaves / (1024.0 * 1024.0),
					AllocatedBytes / (1024.0 * 1024.0));
		for(i = 0; i < NTAB; i++)
		{
			u = 3.0 / NTAB * (i + 0.5);
			shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
			shortrange_table_potential[i] = erfc(u);
#ifdef DISTORTIONTENSORPS
			shortrange_table_tidal[i] = 4.0 * u * u * u / sqrt(M_PI) * exp(-u * u);
#endif
		}
	}
}


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
	if(tree_allocated_flag)
	{
		myfree(Father);
		myfree(Nextnode);
		myfree(Extnodes_base);
		myfree(Nodes_base);
		myfree(DomainNodeIndex);
		tree_allocated_flag = 0;
	}
}




/*! This function does the force computation with direct summation for the
 *  specified particle in the communication buffer. This can be useful for
 *  debugging purposes, in particular for explicit checks of the force
 *  accuracy.
 */
#ifdef FORCETEST
int force_treeevaluate_direct(int target, int mode)
{
	double epsilon, dmax1, dmax2;
	double h, h_inv, dx, dy, dz, r, r2, u, r_inv, fac;
	int i, ptype;
	double pos_x, pos_y, pos_z;
	double acc_x, acc_y, acc_z;

#ifdef PERIODIC
	double fcorr[3];
#endif
#ifdef PERIODIC
	double boxsize, boxhalf;

	boxsize = All.BoxSize;
	boxhalf = 0.5 * All.BoxSize;
#endif
	acc_x = 0;
	acc_y = 0;
	acc_z = 0;
	if(mode == 0)
	{
		pos_x = P[target].Pos[0];
		pos_y = P[target].Pos[1];
		pos_z = P[target].Pos[2];
		ptype = P[target].Type;
	}
	else
	{
		pos_x = GravDataGet[target].Pos[0];
		pos_y = GravDataGet[target].Pos[1];
		pos_z = GravDataGet[target].Pos[2];
#if defined(UNEQUALSOFTENINGS) || defined(SCALARFIELD)
		ptype = GravDataGet[target].Type;
#else
		ptype = P[0].Type;
#endif
	}


	for(i = 0; i < NumPart; i++)
	{
		epsilon = DMAX(All.ForceSoftening[P[i].Type], All.ForceSoftening[ptype]);
		h = epsilon;
		h_inv = 1 / h;
		dx = P[i].Pos[0] - pos_x;
		dy = P[i].Pos[1] - pos_y;
		dz = P[i].Pos[2] - pos_z;
#ifdef PERIODIC
while(dx > boxhalf)
	dx -= boxsize;
while(dy > boxhalf)
	dy -= boxsize;
while(dz > boxhalf)
	dz -= boxsize;
while(dx < -boxhalf)
	dx += boxsize;
while(dy < -boxhalf)
	dy += boxsize;
while(dz < -boxhalf)
	dz += boxsize;
#endif
r2 = dx * dx + dy * dy + dz * dz;
r = sqrt(r2);
u = r * h_inv;
if(u >= 1)
{
	r_inv = 1 / r;
	fac = P[i].Mass * r_inv * r_inv * r_inv;
}
else
{
	if(u < 0.5)
		fac = P[i].Mass * h_inv * h_inv * h_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	else
		fac =
				P[i].Mass * h_inv * h_inv * h_inv * (21.333333333333 -
						48.0 * u + 38.4 * u * u -
						10.666666666667 * u * u *
						u - 0.066666666667 / (u * u * u));
}

acc_x += dx * fac;
acc_y += dy * fac;
acc_z += dz * fac;
#ifdef PERIODIC
if(u > 1.0e-5)
{
	ewald_corr(dx, dy, dz, fcorr);
	acc_x += P[i].Mass * fcorr[0];
	acc_y += P[i].Mass * fcorr[1];
	acc_z += P[i].Mass * fcorr[2];
}
#endif
	}


	if(mode == 0)
	{
		P[target].GravAccelDirect[0] = acc_x;
		P[target].GravAccelDirect[1] = acc_y;
		P[target].GravAccelDirect[2] = acc_z;
	}
	else
	{
		GravDataResult[target].Acc[0] = acc_x;
		GravDataResult[target].Acc[1] = acc_y;
		GravDataResult[target].Acc[2] = acc_z;
	}


	return NumPart;
}
#endif


/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void dump_particles(void)
{
	FILE *fd;
	char buffer[200];
	int i;

	sprintf(buffer, "particles%d.dat", ThisTask);
	fd = fopen(buffer, "w");
	my_fwrite(&NumPart, 1, sizeof(int), fd);
	for(i = 0; i < NumPart; i++)
		my_fwrite(&P[i].Pos[0], 3, sizeof(MyFloat), fd);
	for(i = 0; i < NumPart; i++)
		my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
	for(i = 0; i < NumPart; i++)
		my_fwrite(&P[i].ID, 1, sizeof(int), fd);
	fclose(fd);
}



#ifdef PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin. These corrections are obtained by Ewald summation. (See
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, the Ewald correction is not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization.  The Ewald summation is done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void ewald_init(void)
{
#ifndef NOGRAVITY
	int i, j, k, beg, len, size, n, task, count;
	double x[3], force[3];
	char buf[200];
	FILE *fd;

	if(ThisTask == 0)
	{
		printf("initialize Ewald correction...\n");
		fflush(stdout);
	}

#ifdef DOUBLEPRECISION
	sprintf(buf, "ewald_spc_table_%d_dbl.dat", EN);
#else
	sprintf(buf, "ewald_spc_table_%d.dat", EN);
#endif
	if((fd = fopen(buf, "r")))
	{
		if(ThisTask == 0)
		{
			printf("\nreading Ewald tables from file `%s'\n", buf);
			fflush(stdout);
		}

		my_fread(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
		my_fread(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
		my_fread(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
		my_fread(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
		fclose(fd);
	}
	else
	{
		if(ThisTask == 0)
		{
			printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);
			fflush(stdout);
		}

		/* ok, let's recompute things. Actually, we do that in parallel. */

		size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;
		beg = ThisTask * size;
		len = size;
		if(ThisTask == (NTask - 1))
			len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
		for(i = 0, count = 0; i <= EN; i++)
			for(j = 0; j <= EN; j++)
				for(k = 0; k <= EN; k++)
				{
					n = (i * (EN + 1) + j) * (EN + 1) + k;
					if(n >= beg && n < (beg + len))
					{
						if(ThisTask == 0)
						{
							if((count % (len / 20)) == 0)
							{
								printf("%4.1f percent done\n", count / (len / 100.0));
								fflush(stdout);
							}
						}

						x[0] = 0.5 * ((double) i) / EN;
						x[1] = 0.5 * ((double) j) / EN;
						x[2] = 0.5 * ((double) k) / EN;
						ewald_force(i, j, k, x, force);
						fcorrx[i][j][k] = force[0];
						fcorry[i][j][k] = force[1];
						fcorrz[i][j][k] = force[2];
						if(i + j + k == 0)
							potcorr[i][j][k] = 2.8372975;
						else
							potcorr[i][j][k] = ewald_psi(x);
						count++;
					}
				}

		for(task = 0; task < NTask; task++)
		{
			beg = task * size;
			len = size;
			if(task == (NTask - 1))
				len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
			MPI_Bcast(&fcorrx[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
			MPI_Bcast(&fcorry[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
			MPI_Bcast(&fcorrz[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
			MPI_Bcast(&potcorr[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
		}

		if(ThisTask == 0)
		{
			printf("\nwriting Ewald tables to file `%s'\n", buf);
			fflush(stdout);
			if((fd = fopen(buf, "w")))
			{
				my_fwrite(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
				my_fwrite(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
				my_fwrite(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
				my_fwrite(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
				fclose(fd);
			}
		}
	}

	fac_intp = 2 * EN / All.BoxSize;
	for(i = 0; i <= EN; i++)
		for(j = 0; j <= EN; j++)
			for(k = 0; k <= EN; k++)
			{
				potcorr[i][j][k] /= All.BoxSize;
				fcorrx[i][j][k] /= All.BoxSize * All.BoxSize;
				fcorry[i][j][k] /= All.BoxSize * All.BoxSize;
				fcorrz[i][j][k] /= All.BoxSize * All.BoxSize;
			}

	if(ThisTask == 0)
	{
		printf("initialization of periodic boundaries finished.\n");
		fflush(stdout);
	}
#endif
}


/*! This function looks up the correction force due to the infinite number
 *  of periodic particle/node images. We here use trilinear interpolation
 *  to get it from the precomputed tables, which contain one octant
 *  around the target particle at the origin. The other octants are
 *  obtained from it by exploiting the symmetry properties.
 */
#ifdef FORCETEST
void ewald_corr(double dx, double dy, double dz, double *fper)
{
	int signx, signy, signz;
	int i, j, k;
	double u, v, w;
	double f1, f2, f3, f4, f5, f6, f7, f8;

	if(dx < 0)
	{
		dx = -dx;
		signx = +1;
	}
	else
		signx = -1;
	if(dy < 0)
	{
		dy = -dy;
		signy = +1;
	}
	else
		signy = -1;
	if(dz < 0)
	{
		dz = -dz;
		signz = +1;
	}
	else
		signz = -1;
	u = dx * fac_intp;
	i = (int) u;
	if(i >= EN)
		i = EN - 1;
	u -= i;
	v = dy * fac_intp;
	j = (int) v;
	if(j >= EN)
		j = EN - 1;
	v -= j;
	w = dz * fac_intp;
	k = (int) w;
	if(k >= EN)
		k = EN - 1;
	w -= k;
	f1 = (1 - u) * (1 - v) * (1 - w);
	f2 = (1 - u) * (1 - v) * (w);
	f3 = (1 - u) * (v) * (1 - w);
	f4 = (1 - u) * (v) * (w);
	f5 = (u) * (1 - v) * (1 - w);
	f6 = (u) * (1 - v) * (w);
	f7 = (u) * (v) * (1 - w);
	f8 = (u) * (v) * (w);
	fper[0] = signx * (fcorrx[i][j][k] * f1 +
			fcorrx[i][j][k + 1] * f2 +
			fcorrx[i][j + 1][k] * f3 +
			fcorrx[i][j + 1][k + 1] * f4 +
			fcorrx[i + 1][j][k] * f5 +
			fcorrx[i + 1][j][k + 1] * f6 +
			fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8);
	fper[1] =
			signy * (fcorry[i][j][k] * f1 + fcorry[i][j][k + 1] * f2 +
					fcorry[i][j + 1][k] * f3 + fcorry[i][j + 1][k + 1] * f4 +
					fcorry[i + 1][j][k] * f5 + fcorry[i + 1][j][k + 1] * f6 +
					fcorry[i + 1][j + 1][k] * f7 + fcorry[i + 1][j + 1][k + 1] * f8);
	fper[2] =
			signz * (fcorrz[i][j][k] * f1 + fcorrz[i][j][k + 1] * f2 +
					fcorrz[i][j + 1][k] * f3 + fcorrz[i][j + 1][k + 1] * f4 +
					fcorrz[i + 1][j][k] * f5 + fcorrz[i + 1][j][k + 1] * f6 +
					fcorrz[i + 1][j + 1][k] * f7 + fcorrz[i + 1][j + 1][k + 1] * f8);
}
#endif


/*! This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images. We here use tri-linear
 *  interpolation to get it from the precomputed table, which contains
 *  one octant around the target particle at the origin. The other
 *  octants are obtained from it by exploiting symmetry properties.
 */
double ewald_pot_corr(double dx, double dy, double dz)
{
	int i, j, k;
	double u, v, w;
	double f1, f2, f3, f4, f5, f6, f7, f8;

	if(dx < 0)
		dx = -dx;
	if(dy < 0)
		dy = -dy;
	if(dz < 0)
		dz = -dz;
	u = dx * fac_intp;
	i = (int) u;
	if(i >= EN)
		i = EN - 1;
	u -= i;
	v = dy * fac_intp;
	j = (int) v;
	if(j >= EN)
		j = EN - 1;
	v -= j;
	w = dz * fac_intp;
	k = (int) w;
	if(k >= EN)
		k = EN - 1;
	w -= k;
	f1 = (1 - u) * (1 - v) * (1 - w);
	f2 = (1 - u) * (1 - v) * (w);
	f3 = (1 - u) * (v) * (1 - w);
	f4 = (1 - u) * (v) * (w);
	f5 = (u) * (1 - v) * (1 - w);
	f6 = (u) * (1 - v) * (w);
	f7 = (u) * (v) * (1 - w);
	f8 = (u) * (v) * (w);
	return potcorr[i][j][k] * f1 +
			potcorr[i][j][k + 1] * f2 +
			potcorr[i][j + 1][k] * f3 +
			potcorr[i][j + 1][k + 1] * f4 +
			potcorr[i + 1][j][k] * f5 +
			potcorr[i + 1][j][k + 1] * f6 + potcorr[i + 1][j + 1][k] * f7 + potcorr[i + 1][j + 1][k + 1] * f8;
}



/*! This function computes the potential correction term by means of Ewald
 *  summation.
 */
double ewald_psi(double x[3])
{
	double alpha, psi;
	double r, sum1, sum2, hdotx;
	double dx[3];
	int i, n[3], h[3], h2;

	alpha = 2.0;
	for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
		for(n[1] = -4; n[1] <= 4; n[1]++)
			for(n[2] = -4; n[2] <= 4; n[2]++)
			{
				for(i = 0; i < 3; i++)
					dx[i] = x[i] - n[i];
				r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
				sum1 += erfc(alpha * r) / r;
			}

	for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
		for(h[1] = -4; h[1] <= 4; h[1]++)
			for(h[2] = -4; h[2] <= 4; h[2]++)
			{
				hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
				h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
				if(h2 > 0)
					sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
			}

	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;
	return psi;
}


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
	double alpha, r2;
	double r, val, hdotx, dx[3];
	int i, h[3], n[3], h2;

	alpha = 2.0;
	for(i = 0; i < 3; i++)
		force[i] = 0;
	if(iii == 0 && jjj == 0 && kkk == 0)
		return;
	r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	for(i = 0; i < 3; i++)
		force[i] += x[i] / (r2 * sqrt(r2));
	for(n[0] = -4; n[0] <= 4; n[0]++)
		for(n[1] = -4; n[1] <= 4; n[1]++)
			for(n[2] = -4; n[2] <= 4; n[2]++)
			{
				for(i = 0; i < 3; i++)
					dx[i] = x[i] - n[i];
				r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
				val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);
				for(i = 0; i < 3; i++)
					force[i] -= dx[i] / (r * r * r) * val;
			}

	for(h[0] = -4; h[0] <= 4; h[0]++)
		for(h[1] = -4; h[1] <= 4; h[1]++)
			for(h[2] = -4; h[2] <= 4; h[2]++)
			{
				hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
				h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
				if(h2 > 0)
				{
					val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);
					for(i = 0; i < 3; i++)
						force[i] -= h[i] * val;
				}
			}
}

#endif

#ifdef ADAPTGRAVSOFT
/*! This function updates the gravitational softening of tree nodes.
 *  You need to call it after ags_density() and before the computation of the gravitational force.
 */
void ags_force_update_hmax(void)
{
	int i, no, ta, totDomainNumChanged;
	int *domainList_all;
	int *counts, *offset_list, *offset_hmax;
	MyFloat *domainHmax_loc, *domainHmax_all;

	GlobFlag++;

	DomainNumChanged = 0;
	DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));
#ifndef AGS_UPDATEALLPARTICLES
	for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
#else
		for(i = 0; i < NumPart; i++)
#endif
		{
			no = Father[i];

			while(no >= 0)
			{
				force_drift_node(no, All.Ti_Current);

				if(P[i].AGS_Hsml > Nodes[no].maxsoft)
				{

					Nodes[no].maxsoft = P[i].AGS_Hsml;


					if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node */
					{
						if(Nodes[no].Flag != GlobFlag)
						{
							Nodes[no].Flag = GlobFlag;
							DomainList[DomainNumChanged++] = no;
						}
						break;
					}
				}
				else
					break;

				no = Nodes[no].u.d.father;
			}

		}
	/* share the hmax-data of the pseudo-particles accross CPUs */

	counts = (int *) mymalloc("counts", sizeof(int) * NTask);
	offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
	offset_hmax = (int *) mymalloc("offset_hmax", sizeof(int) * NTask);

	domainHmax_loc = (MyFloat *) mymalloc("domainHmax_loc", DomainNumChanged * 2 * sizeof(MyFloat));

	for(i = 0; i < DomainNumChanged; i++)
	{
		domainHmax_loc[2 * i] = Nodes[DomainList[i]].maxsoft;

	}


	MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

	for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_hmax[0] = 0; ta < NTask; ta++)
	{
		totDomainNumChanged += counts[ta];
		if(ta > 0)
		{
			offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
			offset_hmax[ta] = offset_hmax[ta - 1] + counts[ta - 1] * 2 * sizeof(MyFloat);
		}
	}

	if(ThisTask == 0)
		printf("*ADAPTGRAVSOFT* Hmax exchange: %d topleaves out of %d\n", totDomainNumChanged, NTopleaves);

	domainHmax_all = (MyFloat *) mymalloc("domainHmax_all", totDomainNumChanged * 2 * sizeof(MyFloat));
	domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

	MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
			domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

	for(ta = 0; ta < NTask; ta++)
		counts[ta] *= 2 * sizeof(MyFloat);

	MPI_Allgatherv(domainHmax_loc, 2 * DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
			domainHmax_all, counts, offset_hmax, MPI_BYTE, MPI_COMM_WORLD);


	for(i = 0; i < totDomainNumChanged; i++)
	{
		no = domainList_all[i];

		if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the hmax is updated twice */
			no = Nodes[no].u.d.father;

		while(no >= 0)
		{
			force_drift_node(no, All.Ti_Current);

			if(domainHmax_all[2 * i] > Nodes[no].maxsoft)
			{

				Nodes[no].maxsoft = domainHmax_all[2 * i];

			}
			else
				break;

			no = Nodes[no].u.d.father;
		}
	}


	myfree(domainList_all);
	myfree(domainHmax_all);
	myfree(domainHmax_loc);
	myfree(offset_hmax);
	myfree(offset_list);
	myfree(counts);
	myfree(DomainList);

	CPU_Step[CPU_AGSTREEHMAXUPD] += measure_time();
}
#endif
