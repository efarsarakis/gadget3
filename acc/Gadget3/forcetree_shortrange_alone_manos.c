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
int force_treeevaluate_shortrange(int m_target, int m_mode, int *m_exportflag, int *m_exportnodecount,
		int *m_exportindex)
{
	struct NODE *m_nop = 0;
	int m_no, m_nodesinlist, m_ptype, m_ninteractions, m_nexp, m_tabindex, m_task, m_listindex = 0;
	double m_r2, m_dx, m_dy, m_dz, m_mass, m_r, m_fac, m_u, m_h, m_h_inv, m_h3_inv;
	double m_dxx, m_dyy, m_dzz, m_pdxx, m_pdyy, m_pdzz;
	double m_pos_x, m_pos_y, m_pos_z, m_aold;
	double m_eff_dist;
	double m_rcut, m_asmth, m_asmthfac, m_rcut2, m_dist;
	MyLongDouble m_acc_x, m_acc_y, m_acc_z;
	// cache some global vars in local vars to help compiler with alias analysis
	int m_maxPart = All.MaxPart;
	int m_bunchSize = All.BunchSize;
	int m_maxNodes = MaxNodes;
	integertime m_ti_Current = All.Ti_Current;
	double m_errTol2 = All.ErrTolTheta * All.ErrTolTheta;


	//used/1ST ////////////////////////////////////////////////////////////////////////////////////////////
	double m_xtmp;



	m_acc_x = 0;
	m_acc_y = 0;
	m_acc_z = 0;
	m_ninteractions = 0;
	m_nodesinlist = 0;

	m_rcut = All.Rcut[0];
	m_asmth = All.Asmth[0];

	if(m_mode != 0 && m_mode != 1)
	{
		printf("%d %d %d %d %d\n", m_target, m_mode, *m_exportflag, *m_exportnodecount, *m_exportindex);
		endrun(444);
	}

	if(m_mode == 0)
	{
		m_pos_x = P[m_target].Pos[0];
		m_pos_y = P[m_target].Pos[1];
		m_pos_z = P[m_target].Pos[2];
		m_ptype = P[m_target].Type;
		m_aold = All.ErrTolForceAcc * P[m_target].OldAcc;

	}
	else
	{
		m_pos_x = GravDataGet[m_target].Pos[0];
		m_pos_y = GravDataGet[m_target].Pos[1];
		m_pos_z = GravDataGet[m_target].Pos[2];
		//used/13TH ////////////////////////////////////////////////////////////////////////////////////////////
		m_ptype = P[0].Type;

		m_aold = All.ErrTolForceAcc * GravDataGet[m_target].OldAcc;

	}

	m_rcut2 = m_rcut * m_rcut;

	m_asmthfac = 0.5 / m_asmth * (NTAB / 3.0);

	//used/15TH ////////////////////////////////////////////////////////////////////////////////////////////
	m_h = All.ForceSoftening[m_ptype];
	m_h_inv = 1.0 / m_h;
	m_h3_inv = m_h_inv * m_h_inv * m_h_inv;


	if(m_mode == 0)
	{
		m_no = m_maxPart;		/* root node */
	}
	else
	{
		m_nodesinlist++;
		m_no = GravDataGet[m_target].NodeList[0];
		m_no = Nodes[m_no].u.d.nextnode;	/* open it */
	}

	while(m_no >= 0)
	{
		while(m_no >= 0)
		{
			if(m_no < m_maxPart)
			{
				/* the index of the node is the index of the particle */
				if(P[m_no].Ti_current != m_ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					drift_particle(m_no, m_ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				m_dx = P[m_no].Pos[0] - m_pos_x;
				m_dy = P[m_no].Pos[1] - m_pos_y;
				m_dz = P[m_no].Pos[2] - m_pos_z;

				m_dx = NEAREST(m_dx);
				m_dy = NEAREST(m_dy);
				m_dz = NEAREST(m_dz);

				m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

				m_mass = P[m_no].Mass;

				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;

					P[m_no].GravCost[TakeLevel] += 1.0;
					UNLOCK_WORKCOUNT;
				}
				m_no = Nextnode[m_no];
			}
			else			/* we have an  internal node */
			{
				if(m_no >= m_maxPart + m_maxNodes)	/* pseudo particle */
				{
					if(m_mode == 0)
					{
						if(m_exportflag[m_task = DomainTask[m_no - (m_maxPart + m_maxNodes)]] != m_target)
						{
							m_exportflag[m_task] = m_target;
							m_exportnodecount[m_task] = NODELISTLENGTH;
						}

						if(m_exportnodecount[m_task] == NODELISTLENGTH)
						{
							int m_exitFlag = 0;
							LOCK_NEXPORT;
#pragma omp critical(_nexport_)
							{
								if(Nexport >= m_bunchSize)
								{
									/* out of buffer space. Need to discard work for this particle and interrupt */
									BufferFullFlag = 1;
									m_exitFlag = 1;
								}
								else
								{
									m_nexp = Nexport;
									Nexport++;
								}
							}
							UNLOCK_NEXPORT;
							if(m_exitFlag)
								return -1;

							m_exportnodecount[m_task] = 0;
							m_exportindex[m_task] = m_nexp;
							DataIndexTable[m_nexp].Task = m_task;
							DataIndexTable[m_nexp].Index = m_target;
							DataIndexTable[m_nexp].IndexGet = m_nexp;
						}

						DataNodeList[m_exportindex[m_task]].NodeList[m_exportnodecount[m_task]++] =
								DomainNodeIndex[m_no - (m_maxPart + m_maxNodes)];

						if(m_exportnodecount[m_task] < NODELISTLENGTH)
							DataNodeList[m_exportindex[m_task]].NodeList[m_exportnodecount[m_task]] = -1;
					}
					m_no = Nextnode[m_no - m_maxNodes];
					continue;
				}

				m_nop = &Nodes[m_no];

				if(m_mode == 1)
				{
					if(m_nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
					{
						m_no = -1;
						continue;
					}
				}

				if(!(m_nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
				{
					/* open cell */
					m_no = m_nop->u.d.nextnode;
					continue;
				}

				if(m_nop->Ti_current != m_ti_Current)
				{
					LOCK_PARTNODEDRIFT;
#pragma omp critical(_partnodedrift_)
					force_drift_node(m_no, m_ti_Current);
					UNLOCK_PARTNODEDRIFT;
				}

				m_mass = m_nop->u.d.mass;

				m_dx = m_nop->u.d.s[0] - m_pos_x;
				m_dy = m_nop->u.d.s[1] - m_pos_y;
				m_dz = m_nop->u.d.s[2] - m_pos_z;

				m_dx = NEAREST(m_dx);
				m_dy = NEAREST(m_dy);
				m_dz = NEAREST(m_dz);
				m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;
				//used/21ST ////////////////////////////////////////////////////////////////////////////////////////////
				/* check whether we can stop walking along this branch */
				if(m_r2 > m_rcut2)
				{
					m_eff_dist = m_rcut + 0.5 * m_nop->len;

					m_dist = NEAREST(m_nop->center[0] - m_pos_x);
					if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
					{
						m_no = m_nop->u.d.sibling;
						continue;
					}

					m_dist = NEAREST(m_nop->center[1] - m_pos_y);
					if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
					{
						m_no = m_nop->u.d.sibling;
						continue;
					}

					m_dist = NEAREST(m_nop->center[2] - m_pos_z);
					if(m_dist < -m_eff_dist || m_dist > m_eff_dist)
					{
						m_no = m_nop->u.d.sibling;
						continue;
					}
				}


				if(m_errTol2)	/* check Barnes-Hut opening criterion */
				{
					if(m_nop->len * m_nop->len > m_r2 * m_errTol2)
					{
						/* open cell */
						m_no = m_nop->u.d.nextnode;
						continue;
					}
				}
				else		/* check relative opening criterion */
				{

					//used/23RD ////////////////////////////////////////////////////////////////////////////////////////////
					if(m_mass * m_nop->len * m_nop->len > m_r2 * m_r2 * m_aold)
					{
						/* open cell */
						m_no = m_nop->u.d.nextnode;
						continue;
					}

					/* check in addition whether we lie inside the cell */

					if(fabs(m_nop->center[0] - m_pos_x) < 0.60 * m_nop->len)
					{
						if(fabs(m_nop->center[1] - m_pos_y) < 0.60 * m_nop->len)
						{
							if(fabs(m_nop->center[2] - m_pos_z) < 0.60 * m_nop->len)
							{
								m_no = m_nop->u.d.nextnode;
								continue;
							}
						}
					}

				}


				if(TakeLevel >= 0)
				{
					LOCK_WORKCOUNT;
					m_nop->GravCost += 1.0;
					UNLOCK_WORKCOUNT;
				}

				m_no = m_nop->u.d.sibling;	/* ok, node can be used */

			}

			m_r = sqrt(m_r2);

			//used/31ST ////////////////////////////////////////////////////////////////////////////////////////////
			if(m_r >= m_h)
			{
				m_fac = m_mass / (m_r2 * m_r);
			}

			else
			{

				//used34TH  ////////////////////////////////////////////////////////////////////////////////////////////
				m_u = m_r * m_h_inv;
				if(m_u < 0.5)
					m_fac = m_mass * m_h3_inv * (10.666666666667 + m_u * m_u * (32.0 * m_u - 38.4));
				else
					m_fac =
							m_mass * m_h3_inv * (21.333333333333 - 48.0 * m_u +
									38.4 * m_u * m_u - 10.666666666667 * m_u * m_u * m_u - 0.066666666667 / (m_u * m_u * m_u));


			}

			m_tabindex = (int) (m_asmthfac * m_r);

			if(m_tabindex < NTAB)
			{
				m_fac *= shortrange_table[m_tabindex];

				m_acc_x += FLT(m_dx * m_fac);
				m_acc_y += FLT(m_dy * m_fac);
				m_acc_z += FLT(m_dz * m_fac);
			}

			m_ninteractions++;


		}

		if(m_mode == 1)
		{
			m_listindex++;
			if(m_listindex < NODELISTLENGTH)
			{
				m_no = GravDataGet[m_target].NodeList[m_listindex];
				if(m_no >= 0)
				{
					m_nodesinlist++;
					m_no = Nodes[m_no].u.d.nextnode;	/* open it */
				}
			}
		}
	}


	/* store result at the proper place */
	if(m_mode == 0)
	{
		P[m_target].g.dGravAccel[0] = m_acc_x;
		P[m_target].g.dGravAccel[1] = m_acc_y;
		P[m_target].g.dGravAccel[2] = m_acc_z;
	}
	else
	{
		GravDataResult[m_target].Acc[0] = m_acc_x;
		GravDataResult[m_target].Acc[1] = m_acc_y;
		GravDataResult[m_target].Acc[2] = m_acc_z;
		*m_exportflag = m_nodesinlist;
	}

	return m_ninteractions;
}





