#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SUBFIND
#include "subfind.h"
#include "fof.h"

/* this file processes the local groups in serial mode */


#ifndef MAX_NGB_CHECK
#define MAX_NGB_CHECK 2
#endif

static int *Head, *Next, *Tail, *Len;
static struct cand_dat
{
  int head;
  int len;
  int nsub;
  int rank, subnr, parent;
  int bound_length;
}
 *candidates;




int subfind_process_group_serial(int gr, int Offs)
{
  int i, j, k, p, len, subnr, totlen, ss, ngbs, ndiff, N, head = 0, head_attach, count_cand, len_non_gas;
  int listofdifferent[2], count, prev;
  int ngb_index, part_index, nsubs, rank;
  double SubMass, SubPos[3], SubVel[3], SubCM[3], SubVelDisp, SubVmax, SubVmaxRad, SubSpin[3],
    SubHalfMass, SubMassTab[6];
#ifdef LT_ADD_GAL_TO_SUB
  double SubLum[LT_ADD_GAL_TO_SUB], SubLumAtt[LT_ADD_GAL_TO_SUB], SubLumObs[LT_ADD_GAL_TO_SUB], SubMeanAge,
    SubMeanZStar, SubSFR;
#ifdef DUSTATT
  double Dust[DUSTATT];
#else
  double Dust[1];
#endif
#endif
  MyIDType SubMostBoundID;
  static struct unbind_data *ud;

  while(P[Offs].GrNr != Group[gr].GrNr)
    {
      Offs++;
      if(Offs >= NumPart)
	{
	  printf("don't find a particle for groupnr=%d\n", Group[gr].GrNr);
	  endrun(312);
	}
    }

  N = Group[gr].Len;
  GrNr = Group[gr].GrNr;

  for(i = 0; i < N; i++)
    {
      if(P[Offs + i].GrNr != Group[gr].GrNr)
	{
	  printf
	    ("task=%d, gr=%d: don't have the number of particles for GrNr=%d group-len=%d found=%d before=%d\n",
	     ThisTask, gr, Group[gr].GrNr, N, P[Offs + i].GrNr, P[Offs - 1].GrNr);
	  endrun(312);
	}
    }


  candidates = (struct cand_dat *) mymalloc("candidates", N * sizeof(struct cand_dat));

  Head = (int *) mymalloc("Head", N * sizeof(int));
  Next = (int *) mymalloc("Next", N * sizeof(int));
  Tail = (int *) mymalloc("Tail", N * sizeof(int));
  Len = (int *) mymalloc("Len", N * sizeof(int));
  ud = (struct unbind_data *) mymalloc("ud", N * sizeof(struct unbind_data));

  Head -= Offs;
  Next -= Offs;
  Tail -= Offs;
  Len -= Offs;

  for(i = 0; i < N; i++)
    {
      ud[i].index = Offs + i;
    }

  subfind_loctree_findExtent(N, ud);

  subfind_loctree_treebuild(N, ud);	/* build tree for all particles of this group */

  for(i = Offs; i < Offs + N; i++)
    Head[i] = Next[i] = Tail[i] = -1;

  /* note: particles are already ordered in the order of decreasing density */

  for(i = 0, count_cand = 0; i < N; i++)
    {
      part_index = Offs + i;

      subfind_locngb_treefind(P[part_index].Pos, All.DesLinkNgb, P[part_index].DM_Hsml);

      /* note: returned neighbours are already sorted by distance */

      for(k = 0, ndiff = 0, ngbs = 0; k < All.DesLinkNgb && ngbs < MAX_NGB_CHECK && ndiff < 2; k++)
	{
	  ngb_index = R2list[k].index;

	  if(ngb_index != part_index)	/* to exclude the particle itself */
	    {
	      /* we only look at neighbours that are denser */
	      if(P[ngb_index].u.DM_Density > P[part_index].u.DM_Density)
		{
		  ngbs++;

		  if(Head[ngb_index] >= 0)	/* neighbor is attached to a group */
		    {
		      if(ndiff == 1)
			if(listofdifferent[0] == Head[ngb_index])
			  continue;

		      /* a new group has been found */
		      listofdifferent[ndiff++] = Head[ngb_index];
		    }
		  else
		    {
		      printf("this may not occur.\n");
		      printf
			("ThisTask=%d gr=%d k=%d i=%d part_index=%d ngb_index = %d  head[ngb_index]=%d P[part_index].DM_Density=%g %g GrNrs= %d %d \n",
			 ThisTask, gr, k, i, part_index, ngb_index, Head[ngb_index],
			 P[part_index].u.DM_Density, P[ngb_index].u.DM_Density, P[part_index].GrNr,
			 P[ngb_index].GrNr);
		      endrun(2);
		    }
		}
	    }
	}

      switch (ndiff)		/* treat the different possible cases */
	{
	case 0:		/* this appears to be a lonely maximum -> new group */
	  head = part_index;
	  Head[part_index] = Tail[part_index] = part_index;
	  Len[part_index] = 1;
	  Next[part_index] = -1;
	  break;

	case 1:		/* the particle is attached to exactly one group */
	  head = listofdifferent[0];
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;
	  break;

	case 2:		/* the particle merges two groups together */

	  head = listofdifferent[0];
	  head_attach = listofdifferent[1];

	  if(Len[head_attach] > Len[head])	/* other group is longer, swap them */
	    {
	      head = listofdifferent[1];
	      head_attach = listofdifferent[0];
	    }

	  /* only in case the attached group is long enough we bother to register is 
	     as a subhalo candidate */

	  if(Len[head_attach] >= All.DesLinkNgb)
	    {
	      candidates[count_cand].len = Len[head_attach];
	      candidates[count_cand].head = Head[head_attach];
	      count_cand++;
	    }

	  /* now join the two groups */
	  Next[Tail[head]] = head_attach;
	  Tail[head] = Tail[head_attach];
	  Len[head] += Len[head_attach];

	  ss = head_attach;
	  do
	    {
	      Head[ss] = head;
	    }
	  while((ss = Next[ss]) >= 0);

	  /* finally, attach the particle */
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;
	  break;

	default:
	  printf("can't be! (a)\n");
	  endrun(1);
	  break;
	}
    }

  /* add the full thing as a subhalo candidate */
  for(i = 0, prev = -1; i < N; i++)
    {
      if(Head[Offs + i] == Offs + i)
	if(Next[Tail[Offs + i]] == -1)
	  {
	    if(prev < 0)
	      head = Offs + i;
	    if(prev >= 0)
	      Next[prev] = Offs + i;

	    prev = Tail[Offs + i];
	  }
    }

  candidates[count_cand].len = N;
  candidates[count_cand].head = head;
  count_cand++;

  /* go through them once and assign the rank */
  for(i = 0, p = head, rank = 0; i < N; i++)
    {
      Len[p] = rank++;
      p = Next[p];
    }

  /* for each candidate, we now pull out the rank of its head */
  for(k = 0; k < count_cand; k++)
    candidates[k].rank = Len[candidates[k].head];

  for(i = Offs; i < Offs + N; i++)
    Tail[i] = -1;

  for(k = 0, nsubs = 0; k < count_cand; k++)
    {
      for(i = 0, p = candidates[k].head, len = 0; i < candidates[k].len; i++, p = Next[p])
	if(Tail[p] < 0)
	  ud[len++].index = p;

      if(len >= All.DesLinkNgb)
	len = subfind_unbind(ud, len, &len_non_gas);

#ifdef NO_GAS_CLOUDS
      if(len_non_gas >= All.DesLinkNgb)
#else
      if(len >= All.DesLinkNgb)
#endif
	{
	  /* ok, we found a substructure */

	  for(i = 0; i < len; i++)
	    Tail[ud[i].index] = nsubs;	/* we use this to flag the substructures */

	  candidates[k].nsub = nsubs;
	  candidates[k].bound_length = len;
	  nsubs++;
	}
      else
	{
	  candidates[k].nsub = -1;
	  candidates[k].bound_length = 0;
	}
    }

#ifdef VERBOSE
  printf("\nGroupLen=%d  (gr=%d)\n", N, gr);
  printf("Number of substructures: %d\n", nsubs);
#endif

  Group[gr].Nsubs = nsubs;
  Group[gr].Pos[0] = Group[gr].CM[0];
  Group[gr].Pos[1] = Group[gr].CM[1];
  Group[gr].Pos[2] = Group[gr].CM[2];

#ifdef OMP_SORT
  omp_qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_boundlength);
#else
  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_boundlength);
#endif

  /* now we determine the parent subhalo for each candidate */
  for(k = 0; k < count_cand; k++)
    {
      candidates[k].subnr = k;
      candidates[k].parent = 0;
    }

#ifdef OMP_SORT
  omp_qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_rank);
#else
  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_rank);
#endif

  for(k = 0; k < count_cand; k++)
    {
      for(j = k + 1; j < count_cand; j++)
	{
	  if(candidates[j].rank > candidates[k].rank + candidates[k].len)
	    break;

	  if(candidates[k].rank + candidates[k].len >= candidates[j].rank + candidates[j].len)
	    {
	      if(candidates[k].bound_length >= All.DesLinkNgb)
		candidates[j].parent = candidates[k].subnr;
	    }
	  else
	    {
	      printf("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
		     k, count_cand, (int) candidates[k].rank, candidates[k].len,
		     (int) candidates[k].bound_length, candidates[j].rank,
		     (int) candidates[j].len, candidates[j].bound_length);
	      endrun(121235513);
	    }
	}
    }

#ifdef OMP_SORT
  omp_qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_subnr);
#else
  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_subnr);
#endif
  /* now determine the properties */

  for(k = 0, subnr = 0, totlen = 0; k < nsubs; k++)
    {
      len = candidates[k].bound_length;

#ifdef VERBOSE
      printf("subnr=%d  SubLen=%d\n", subnr, len);
#endif

      totlen += len;

      for(i = 0, p = candidates[k].head, count = 0; i < candidates[k].len; i++)
	{
	  if(Tail[p] == candidates[k].nsub)
	    ud[count++].index = p;

	  p = Next[p];
	}

      if(count != len)
	endrun(12);


#ifdef LT_ADD_GAL_TO_SUB
      subfind_determine_sub_halo_properties(ud, len, &SubMass,
					    &SubPos[0], &SubVel[0], &SubCM[0], &SubVelDisp, &SubVmax,
					    &SubVmaxRad, &SubSpin[0], &SubMostBoundID, &SubHalfMass,
					    &SubMassTab[0], &SubLum[0], &SubLumAtt[0], &SubLumObs[0],
					    &Dust[0], &SubMeanAge, &SubMeanZStar, &SubSFR);
#else
      subfind_determine_sub_halo_properties(ud, len, &SubMass,
					    &SubPos[0], &SubVel[0], &SubCM[0], &SubVelDisp, &SubVmax,
					    &SubVmaxRad, &SubSpin[0], &SubMostBoundID, &SubHalfMass,
					    &SubMassTab[0]);
#endif

      if(Nsubgroups >= MaxNsubgroups)
	endrun(899);

      if(subnr == 0)
	{
	  for(j = 0; j < 3; j++)
	    Group[gr].Pos[j] = SubPos[j];
	}

      SubGroup[Nsubgroups].Len = len;
      if(subnr == 0)
	SubGroup[Nsubgroups].Offset = Group[gr].Offset;
      else
	SubGroup[Nsubgroups].Offset = SubGroup[Nsubgroups - 1].Offset + SubGroup[Nsubgroups - 1].Len;
      SubGroup[Nsubgroups].GrNr = GrNr - 1;
      SubGroup[Nsubgroups].SubNr = subnr;
      SubGroup[Nsubgroups].SubParent = candidates[k].parent;
      SubGroup[Nsubgroups].Mass = SubMass;
      SubGroup[Nsubgroups].SubMostBoundID = SubMostBoundID;
      SubGroup[Nsubgroups].SubVelDisp = SubVelDisp;
      SubGroup[Nsubgroups].SubVmax = SubVmax;
      SubGroup[Nsubgroups].SubVmaxRad = SubVmaxRad;
      SubGroup[Nsubgroups].SubHalfMass = SubHalfMass;

      for(j = 0; j < 3; j++)
	{
	  SubGroup[Nsubgroups].Pos[j] = SubPos[j];
	  SubGroup[Nsubgroups].CM[j] = SubCM[j];
	  SubGroup[Nsubgroups].Vel[j] = SubVel[j];
	  SubGroup[Nsubgroups].Spin[j] = SubSpin[j];
	}

#ifdef SAVE_MASS_TAB
      for(j = 0; j < 6; j++)
	SubGroup[Nsubgroups].MassTab[j] = SubMassTab[j];
#endif

#ifdef LT_ADD_GAL_TO_SUB
      if(SubMassTab[4] > 0.)
	{
	  for(j = 0; j < LT_ADD_GAL_TO_SUB; j++)
	    {
	      if(SubLum[j] > 0)
		SubGroup[Nsubgroups].SubLum[j] = -2.5 * log10(SubLum[j]);
	      else
		SubGroup[Nsubgroups].SubLum[j] = 0;
#ifdef DUSTATT
	      if(SubLumAtt[j] > 0)
		SubGroup[Nsubgroups].SubLumAtt[j] = -2.5 * log10(SubLumAtt[j]);
	      else
		SubGroup[Nsubgroups].SubLumAtt[j] = 0;
#endif
#ifdef OBSERVER_FRAME
	      if(SubLumObs[j] > 0)
		SubGroup[Nsubgroups].SubLumObs[j] = -2.5 * log10(SubLumObs[j]);
	      else
		SubGroup[Nsubgroups].SubLumObs[j] = 0;
#endif
	    }
#ifdef DUSTATT
	  for(j = 0; j < DUSTATT; j++)
	    SubGroup[Nsubgroups].DustProfile[j] = Dust[j];
#endif
	  SubGroup[Nsubgroups].SubMeanAge = SubMeanAge;
	  SubGroup[Nsubgroups].SubMeanZStar = SubMeanZStar;
	}
      else
	{
	  for(j = 0; j < LT_ADD_GAL_TO_SUB; j++)
	    {
	      SubGroup[Nsubgroups].SubLum[j] = 0;
#ifdef DUSTATT
	      SubGroup[Nsubgroups].SubLumAtt[j] = 0;
#endif
#ifdef OBSERVER_FRAME
	      SubGroup[Nsubgroups].SubLumObs[j] = 0;
#endif
	    }
	  SubGroup[Nsubgroups].SubMeanAge = 0;
	  SubGroup[Nsubgroups].SubMeanZStar = 0;
#ifdef DUSTATT
	  for(j = 0; j < DUSTATT; j++)
	    SubGroup[Nsubgroups].DustProfile[j] = 0;
#endif
	}

      if(SubMassTab[0] > 0.)
	{
	  SubGroup[Nsubgroups].SubSFR = SubSFR;
	}
      else
	{
	  SubGroup[Nsubgroups].SubSFR = 0;;
	}
#endif

      Nsubgroups++;

      /* Let's now assign the subgroup number */

      for(i = 0; i < len; i++)
	P[ud[i].index].SubNr = subnr;

      subnr++;
    }

#ifdef VERBOSE
  printf("Fuzz=%d\n", N - totlen);
#endif

  myfree(ud);
  myfree(Len + Offs);
  myfree(Tail + Offs);
  myfree(Next + Offs);
  myfree(Head + Offs);

  myfree(candidates);

  return Offs;
}




int subfind_unbind(struct unbind_data *ud, int len, int *len_non_gas)
{
  double *bnd_energy, energy_limit, weakly_bound_limit = 0;
  int i, j, p, minindex, unbound, phaseflag;
  double ddxx, s[3], dx[3], v[3], dv[3], pos[3];
  double vel_to_phys, H_of_a, atime, pot, minpot = 0;
  double boxsize;
  double TotMass;

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      vel_to_phys = atime = 1;
      H_of_a = 0;
    }

  bnd_energy = (double *) mymalloc("bnd_energy", len * sizeof(double));

  phaseflag = 0;		/* this means we will recompute the potential for all particles */

  do
    {
      subfind_loctree_treebuild(len, ud);

      /* let's compute the potential  */

      if(phaseflag == 0)	/* redo it for all the particles */
	{
	  for(i = 0, minindex = -1, minpot = 1.0e30; i < len; i++)
	    {
	      p = ud[i].index;

	      pot = subfind_loctree_treeevaluate_potential(p);
	      /* note: add self-energy */
#ifndef ADAPTGRAVSOFT
	      P[p].u.DM_Potential = pot + P[p].Mass / All.SofteningTable[P[p].Type];
#else
	      P[p].u.DM_Potential = pot + P[p].Mass * 2.8 / P[p].AGS_Hsml;
#endif
	      P[p].u.DM_Potential *= All.G / atime;

	      if(All.TotN_gas > 0 && (FOF_PRIMARY_LINK_TYPES & 1) == 0 &&
		 (FOF_SECONDARY_LINK_TYPES & 1) == 0 && All.OmegaBaryon > 0)
		P[p].u.DM_Potential *= All.Omega0 / (All.Omega0 - All.OmegaBaryon);

	      if(P[p].u.DM_Potential < minpot || minindex == -1)
		{
		  minpot = P[p].u.DM_Potential;
		  minindex = p;
		}
	    }

	  for(j = 0; j < 3; j++)
	    pos[j] = P[minindex].Pos[j];	/* position of minimum potential */
	}
      else
	{
	  /* we only repeat for those close to the unbinding threshold */
	  for(i = 0; i < len; i++)
	    {
	      p = ud[i].index;

	      if(P[p].v.DM_BindingEnergy >= weakly_bound_limit)
		{
		  pot = subfind_loctree_treeevaluate_potential(p);
		  /* note: add self-energy */
#ifndef ADAPTGRAVSOFT
		  P[p].u.DM_Potential = pot + P[p].Mass / All.SofteningTable[P[p].Type];
#else
		  P[p].u.DM_Potential = pot + P[p].Mass * 2.8 / P[p].AGS_Hsml;
#endif
		  P[p].u.DM_Potential *= All.G / atime;

		  if(All.TotN_gas > 0 && (FOF_PRIMARY_LINK_TYPES & 1) == 0 &&
		     (FOF_SECONDARY_LINK_TYPES & 1) == 0 && All.OmegaBaryon > 0)
		    P[p].u.DM_Potential *= All.Omega0 / (All.Omega0 - All.OmegaBaryon);
		}
	    }
	}

      /* let's get bulk velocity and the center-of-mass */

      v[0] = v[1] = v[2] = 0;
      s[0] = s[1] = s[2] = 0;

      for(i = 0, TotMass = 0; i < len; i++)
	{
	  p = ud[i].index;

	  for(j = 0; j < 3; j++)
	    {
#ifdef PERIODIC
	      ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	      ddxx = P[p].Pos[j] - pos[j];
#endif
	      s[j] += P[p].Mass * ddxx;
	      v[j] += P[p].Mass * P[p].Vel[j];
	    }
	  TotMass += P[p].Mass;
	}

      for(j = 0; j < 3; j++)
	{
	  v[j] /= TotMass;
	  s[j] /= TotMass;	/* center-of-mass */

	  s[j] += pos[j];

#ifdef PERIODIC
	  while(s[j] < 0)
	    s[j] += boxsize;
	  while(s[j] >= boxsize)
	    s[j] -= boxsize;
#endif
	}

      for(i = 0; i < len; i++)
	{
	  p = ud[i].index;

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
#ifdef PERIODIC
	      dx[j] = atime * NEAREST(P[p].Pos[j] - s[j]);
#else
	      dx[j] = atime * (P[p].Pos[j] - s[j]);
#endif
	      dv[j] += H_of_a * dx[j];
	    }

	  P[p].v.DM_BindingEnergy =
	    P[p].u.DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

#ifdef DENSITY_SPLIT_BY_TYPE
	  if(P[p].Type == 0)
	    P[p].v.DM_BindingEnergy += P[p].w.int_energy;
#endif
	  bnd_energy[i] = P[p].v.DM_BindingEnergy;
	}

#ifdef OMP_SORT
      omp_qsort(bnd_energy, len, sizeof(double), subfind_compare_binding_energy);	/* largest comes first! */
#else
      qsort(bnd_energy, len, sizeof(double), subfind_compare_binding_energy);	/* largest comes first! */
#endif

      energy_limit = bnd_energy[(int) (0.25 * len)];

      for(i = 0, unbound = 0; i < len - 1; i++)
	{
	  if(bnd_energy[i] > 0)
	    unbound++;
	  else
	    unbound--;

	  if(unbound <= 0)
	    break;
	}
      weakly_bound_limit = bnd_energy[i];

      /* now omit unbound particles,  but at most 1/4 of the original size */

      for(i = 0, unbound = 0, *len_non_gas = 0; i < len; i++)
	{
	  p = ud[i].index;
	  if(P[p].v.DM_BindingEnergy > 0 && P[p].v.DM_BindingEnergy > energy_limit)
	    {
	      unbound++;
	      ud[i] = ud[len - 1];
	      i--;
	      len--;
	    }
	  else if(P[p].Type != 0)
	    (*len_non_gas)++;
	}

      if(len < All.DesLinkNgb)
	break;

      if(phaseflag == 0)
	{
	  if(unbound > 0)
	    phaseflag = 1;
	}
      else
	{
	  if(unbound == 0)
	    {
	      phaseflag = 0;	/* this will make us repeat everything once more for all particles */
	      unbound = 1;
	    }
	}
    }
  while(unbound > 0);

  myfree(bnd_energy);

  return (len);
}



int subfind_compare_grp_particles(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < ((struct particle_data *) b)->GrNr)
    return -1;

  if(((struct particle_data *) a)->GrNr > ((struct particle_data *) b)->GrNr)
    return +1;

  if(((struct particle_data *) a)->SubNr < ((struct particle_data *) b)->SubNr)
    return -1;

  if(((struct particle_data *) a)->SubNr > ((struct particle_data *) b)->SubNr)
    return +1;

  if(((struct particle_data *) a)->v.DM_BindingEnergy < ((struct particle_data *) b)->v.DM_BindingEnergy)
    return -1;

  if(((struct particle_data *) a)->v.DM_BindingEnergy > ((struct particle_data *) b)->v.DM_BindingEnergy)
    return +1;

  return 0;
}

#ifdef LT_ADD_GAL_TO_SUB

void load_CB07_table(int IMFt, int num)
{
  /* This routine reads the CB07 tables for luminosities */

  FILE *TempiFile, *CB07File, *Filters_Effective_WavelenghtFile;
  int imag = 0, imet = 0, i;
  char buff[200];

  sprintf(buff, "%s/Filters_Effective_Wavelenght.txt", All.BC_SED_Path);
  if(!(Filters_Effective_WavelenghtFile = fopen(buff, "rb")))
    {
      printf("Could not open file %s' !\n", buff);
      endrun(9540);
    }

  for(i = 0; i < LT_ADD_GAL_TO_SUB; i++)
    {
      fscanf(Filters_Effective_WavelenghtFile, "%f", &Filters_Effective_Wavelenght[i]);
    }

  fclose(Filters_Effective_WavelenghtFile);

  sprintf(buff, "%s/tempi", All.BC_SED_Path);
  if(!(TempiFile = fopen(buff, "rb")))
    {
      printf("Could not open file %s' !\n", buff);
      endrun(9540);
    }

  for(imet = 0; imet <= 6; imet++)
    {
      sprintf(buff, "%s/%s/cb2007_hr_stelib_m%d2_ssp.multi_mag_vega", All.BC_SED_Path, SFs[0].IMFname,
	      2 + imet);
      if(!(CB07File = fopen(buff, "rb")))
	{
	  printf("Could not open file %s' !\n", buff);
	  endrun(9541);

	}

      for(i = 0; i < 219; i++)
	{
	  if(imet == 0)
	    {
	      fscanf(TempiFile, "%f", &tempiAS[i]);
	      tempiAS[i] = pow(10, tempiAS[i]) / 1.e9;
	    }
	  /* now tempiAS should be in Gyrs */

	  for(imag = 0; imag < LT_ADD_GAL_TO_SUB; ++imag)
	    {
	      fscanf(CB07File, "%f",
		     &CB07[imag + (i * LT_ADD_GAL_TO_SUB) + (imet * 219 * LT_ADD_GAL_TO_SUB)]);
	    }
	}
      fclose(CB07File);

#ifdef OBSERVER_FRAME
      sprintf(buff, "%s/%s/cb2007_hr_stelib_m%d2_ssp.multi_mag_vega_%03d", All.BC_SED_Path, SFs[0].IMFname,
	      2 + imet, num);
      if(!(CB07File = fopen(buff, "rb")))
	{
	  printf("Could not open file %s' !\n", buff);
	  endrun(9542);

	}

      for(i = 0; i < 219; i++)
	for(imag = 0; imag < LT_ADD_GAL_TO_SUB; ++imag)
	  fscanf(CB07File, "%f", &CB07obs[imag + (i * LT_ADD_GAL_TO_SUB) + (imet * 219 * LT_ADD_GAL_TO_SUB)]);

      fclose(CB07File);
#endif
    }

  fclose(TempiFile);
}

void add_luminosity(long PART_index, int *FLAG_AGE, double *lum, double *lum_obs)
{
  float metCB07[8];

  int ind_T, ind_Z, j, k, i;
  float rappt1, rappt, rappz1, rappz, Tstar, Tdiff, Zstar;

  float Lumzt, Lumz1t, Lumzt1, Lumz1t1, Tnow;
#ifdef OBSERVER_FRAME
  float Lumzto, Lumz1to, Lumzt1o, Lumz1t1o;
#endif

  metCB07[0] = 0.0001 / 0.02;	// Normalized to Z_sun = 0.02
  metCB07[1] = 0.0004 / 0.02;
  metCB07[2] = 0.004 / 0.02;
  metCB07[3] = 0.008 / 0.02;
  metCB07[4] = 0.02 / 0.02;
  metCB07[5] = 0.0500 / 0.02;
  metCB07[6] = 0.1000 / 0.02;
  metCB07[7] = 1.e9;

  Tnow = get_age(All.Time);

  Tstar = get_age(P[PART_index].StellarAge);
  Tdiff = Tstar - Tnow;
  if(Tdiff > 1.e-2)
    *FLAG_AGE = 0;
  else
    *FLAG_AGE = 1;


  Zstar = get_metallicity_subfind(PART_index) / 0.02;
  j = 0;
  while(Tdiff >= tempiAS[j])
    j++;
  ind_T = j;
  if(j > 0)
    {
      rappt1 = (tempiAS[ind_T] - Tdiff) / (tempiAS[ind_T] - tempiAS[ind_T - 1]);
      rappt = (Tdiff - tempiAS[ind_T - 1]) / (tempiAS[ind_T] - tempiAS[ind_T - 1]);
    }
  if(j == 0)
    {
      rappt1 = 0.;
      rappt = 1.;
      ind_T++;
    }
  k = 0;
  while(Zstar >= metCB07[k])
    k++;
  ind_Z = k;
  if(k > 0)
    {
      rappz1 = (metCB07[ind_Z] - Zstar) / (metCB07[ind_Z] - metCB07[ind_Z - 1]);
      rappz = (Zstar - metCB07[ind_Z - 1]) / (metCB07[ind_Z] - metCB07[ind_Z - 1]);
    }
  else
    {
      rappz = 1.;
      rappz1 = 0.;
      ind_Z++;
    }


  for(i = 0; i < LT_ADD_GAL_TO_SUB; i++)
    {
      if(ind_Z == 7)
	{
	  rappz1 = 0.;
	  rappz = 1.;
	  Lumz1t = 0.;
	  Lumz1t1 = 0.;
	}
      else
	{
	  Lumz1t = -CB07[(ind_Z) * 219 * LT_ADD_GAL_TO_SUB + (ind_T - 1) * LT_ADD_GAL_TO_SUB + i];
	  Lumz1t1 = -CB07[(ind_Z) * 219 * LT_ADD_GAL_TO_SUB + (ind_T) * LT_ADD_GAL_TO_SUB + i];
#ifdef OBSERVER_FRAME
	  Lumz1to = -CB07obs[(ind_Z) * 219 * LT_ADD_GAL_TO_SUB + (ind_T - 1) * LT_ADD_GAL_TO_SUB + i];
	  Lumz1t1o = -CB07obs[(ind_Z) * 219 * LT_ADD_GAL_TO_SUB + (ind_T) * LT_ADD_GAL_TO_SUB + i];
#endif
	}
      Lumzt = -CB07[(ind_Z - 1) * 219 * LT_ADD_GAL_TO_SUB + (ind_T - 1) * LT_ADD_GAL_TO_SUB + i];
      Lumzt1 = -CB07[(ind_Z - 1) * 219 * LT_ADD_GAL_TO_SUB + (ind_T) * LT_ADD_GAL_TO_SUB + i];
#ifdef OBSERVER_FRAME
      Lumzto = -CB07[(ind_Z - 1) * 219 * LT_ADD_GAL_TO_SUB + (ind_T - 1) * LT_ADD_GAL_TO_SUB + i];
      Lumzt1o = -CB07[(ind_Z - 1) * 219 * LT_ADD_GAL_TO_SUB + (ind_T) * LT_ADD_GAL_TO_SUB + i];
#endif
      lum[i] = ((pow(10, Lumzt / 2.5) * rappz + pow(10, Lumz1t / 2.5) * rappz1) * rappt +
		(pow(10, Lumzt1 / 2.5) * rappz + pow(10, Lumz1t1 / 2.5) * rappz1) * rappt1)
	* MetP[P[PART_index].pt.MetID].iMass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
#ifdef OBSERVER_FRAME
      lum_obs[i] = ((pow(10, Lumzto / 2.5) * rappz + pow(10, Lumz1to / 2.5) * rappz1) * rappt +
		    (pow(10, Lumzt1o / 2.5) * rappz + pow(10, Lumz1t1o / 2.5) * rappz1) * rappt1)
	* MetP[P[PART_index].pt.MetID].iMass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
#endif
    }

}

#endif

#ifdef LT_ADD_GAL_TO_SUB
void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					   double *pos, double *vel, double *cm, double *veldisp,
					   double *vmax, double *vmaxrad, double *spin,
					   MyIDType * mostboundid, double *halfmassrad, double *mass_tab,
					   double *lum, double *lum_att, double *lum_obs, double *dust,
					   double *meanage, double *meanzstar, double *SFRgas)
#else
void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					   double *pos, double *vel, double *cm, double *veldisp,
					   double *vmax, double *vmaxrad, double *spin,
					   MyIDType * mostboundid, double *halfmassrad, double *mass_tab)
#endif
{
  int i, j, p, num_use, i_use;
  double s[3], v[3], max, vel_to_phys, H_of_a, atime, minpot;
  double lx, ly, lz, dv[3], dx[3], disp, rr_tmp, disp_tmp;
  double boxsize, ddxx;
  sort_r2list *rr_list = 0;
  int minindex;
  double mass, maxrad;
  int nstar = 0, ndm = 0, ngas = 0, nbh = 0;

#ifdef LT_ADD_GAL_TO_SUB
  float Temperature, xclouds;
  sort_r2list *rr_list_gal = 0;
  int istar = 0, k;
  int FLAG_AGE;
  double lum_particle[LT_ADD_GAL_TO_SUB], r_lum = 0, lum_particle_obs[LT_ADD_GAL_TO_SUB];
  int ilum;
  *meanage = 0;
  *meanzstar = 0;
  *SFRgas = 0;
  for(k = 0; k < LT_ADD_GAL_TO_SUB; k++)
    {
      lum[k] = 0;
#ifdef OBSERVER_FRAME
      lum_obs[k] = 0;
#endif
    }

#ifdef DUSTATT
  int rad_bin, rad_bin_z;
  double mcold, Mcold_in_bin[DUSTATT], Zgas_in_bin[DUSTATT], TauV_in_bin[DUSTATT];
  double mu, Unit_corrections_for_dust, Tau_Lambda_over_V, rr2_tmp, depth;
#ifdef OBSERVER_FRAME
  double Tau_Lambda_over_V_obs;
#endif
  long mu_seed;
  float MUWIDTH = 0.2, MUCENTER = 0.3;
  float term_depth1 = 0, term_depth2 = 0, term_depth_total = 0;

  for(i = 0; i < LT_ADD_GAL_TO_SUB; i++)
    lum_att[i] = 0;

  for(i = 0; i < DUSTATT; i++)
    {
      Mcold_in_bin[i] = 0;
      Zgas_in_bin[i] = 0;
      TauV_in_bin[i] = 0;
    }

  Unit_corrections_for_dust =
    All.UnitMass_in_g / All.HubbleParam / (1.4 * PROTONMASS * (All.UnitLength_in_cm / All.HubbleParam) *
					   (All.UnitLength_in_cm / All.HubbleParam) * 2.1 * 1.e21);

  mu_seed = -150;

#endif
#endif

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      vel_to_phys = atime = 1;
      H_of_a = 0;
    }

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(P[p].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[p].u.DM_Potential;
	  minindex = p;
	}
      switch (P[p].Type)
	{
	case 0:
	  ngas++;
	  break;
	case 1:
	case 2:
	case 3:
	  ndm++;
	  break;
	case 4:
	  nstar++;
	  break;
	case 5:
	  nbh++;
	  break;
	}
    }

  if(minindex == -1)
    endrun(875412);

  for(j = 0; j < 3; j++)
    pos[j] = P[minindex].Pos[j];


  /* pos[] now holds the position of minimum potential */
  /* we take it that as the center */


  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
#ifdef NO_GAS_CLOUDS
      if(P[p].Type > 0)
#endif
	if(P[p].v.DM_BindingEnergy < minpot || minindex == -1)
	  {
	    minpot = P[p].v.DM_BindingEnergy;
	    minindex = p;
	  }
    }

  if(minindex == -1)
    endrun(875413);

  *mostboundid = P[minindex].ID;


  /* let's get bulk velocity and the center-of-mass */
  /* here we still can take all particles */

  for(j = 0; j < 3; j++)
    s[j] = v[j] = 0;

  for(j = 0; j < 6; j++)
    mass_tab[j] = 0;

  for(i = 0, mass = 0; i < num; i++)
    {
      p = d[i].index;
      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	  ddxx = P[p].Pos[j] - pos[j];
#endif
	  s[j] += P[p].Mass * ddxx;
	  v[j] += P[p].Mass * P[p].Vel[j];
	}
      mass += P[p].Mass;

      mass_tab[P[p].Type] += P[p].Mass;
    }

  *totmass = mass;

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] += pos[j];

#ifdef PERIODIC
      while(s[j] < 0)
	s[j] += boxsize;
      while(s[j] >= boxsize)
	s[j] -= boxsize;
#endif

      cm[j] = s[j];
    }


  disp = lx = ly = lz = 0;

  /* Here we have to perform only on the dm particles for consistency */
  mass = 0;
  num_use = num;
#ifdef DENSITY_SPLIT_BY_TYPE
  num_use = ndm;
#endif

#ifdef LT_ADD_GAL_TO_SUB
  rr_list_gal = (sort_r2list *) mymalloc("rr_list_gal", sizeof(sort_r2list) * num);
#endif
  if(num_use > 0)
    rr_list = (sort_r2list *) mymalloc("rr_list", sizeof(sort_r2list) * num_use);

  for(i = 0, i_use = 0; i < num; i++)
    {
      p = d[i].index;

#ifdef LT_ADD_GAL_TO_SUB
      rr_list_gal[i].index = p;
      rr_list_gal[i].mass = P[p].Mass;
#ifdef DUSTATT
      rr2_tmp = 0;
#endif
#endif

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - s[j]);
#else
	  ddxx = P[p].Pos[j] - s[j];
#endif
	  dx[j] = atime * ddxx;
	  dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];

	  disp_tmp += P[p].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
#ifdef PERIODIC
	  ddxx = NEAREST(P[p].Pos[j] - pos[j]);
#else
	  ddxx = P[p].Pos[j] - pos[j];
#endif
	  ddxx = atime * ddxx;
	  rr_tmp += ddxx * ddxx;
#ifdef LT_ADD_GAL_TO_SUB
#ifdef DUSTATT
	  if(j < 2)
	    rr2_tmp += ddxx * ddxx;
	  if(j == 2)
	    rr_list_gal[i].z = ddxx;
#endif
#endif

	}

      rr_tmp = sqrt(rr_tmp);
#ifdef LT_ADD_GAL_TO_SUB
      rr_list_gal[i].r = rr_tmp;
#ifdef DUSTATT
      rr_list_gal[i].r2 = sqrt(rr2_tmp);
#endif
#endif

#ifdef DENSITY_SPLIT_BY_TYPE
      if(P[p].Type >= 1 && P[p].Type <= 3)  /*-- only for dm part --*/
	{
#endif
	  rr_list[i_use].mass = P[p].Mass;
	  rr_list[i_use].r = rr_tmp;
	  disp += disp_tmp;
	  mass += P[p].Mass;
	  lx += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
	  ly += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
	  lz += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

	  i_use++;
#ifdef DENSITY_SPLIT_BY_TYPE
	}
#endif
    }

  if(i_use != num_use)
    endrun(564321);

  if(num_use > 0)
    {
      *veldisp = sqrt(disp / (3 * mass));	/* convert to 1d velocity dispersion */

      spin[0] = lx / mass;
      spin[1] = ly / mass;
      spin[2] = lz / mass;

#ifdef OMP_SORT
      omp_qsort(rr_list, num_use, sizeof(sort_r2list), subfind_compare_dist_rotcurve);
#else
      qsort(rr_list, num_use, sizeof(sort_r2list), subfind_compare_dist_rotcurve);
#endif

/*--- Here we still have to fix for individual masses, 
        maybe we even want the total mass for this ? 
        Note however that even within multi mass DM simulation
        the DM mass within a clean halo should be the same ... ------*/
      *halfmassrad = rr_list[num_use / 2].r;

      /* compute cumulative mass */
      for(i = 1; i < num_use; i++)
	rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

/*--- Note that here we might want to correct for the baryon fraction ?? ---*/
      for(i = num_use - 1, max = 0, maxrad = 0; i > 5; i--)
	if(rr_list[i].mass / rr_list[i].r > max)
	  {
	    max = rr_list[i].mass / rr_list[i].r;
	    maxrad = rr_list[i].r;
	  }

      *vmax = sqrt(All.G * max);
      *vmaxrad = maxrad;

      myfree(rr_list);
    }
  else
    *veldisp = *halfmassrad = *vmax = *vmaxrad = spin[0] = spin[1] = spin[2] = 0;


#ifdef LT_ADD_GAL_TO_SUB
  istar = 0;

#ifdef OMP_SORT
  omp_qsort(rr_list_gal, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);
#else
  qsort(rr_list_gal, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);
#endif

  /* compute luminosity radius */
  for(i = 0; i < num; i++)
    {
      if(P[rr_list_gal[i].index].Type == 4)
	{
	  istar++;
	  if(nstar == 1)
	    r_lum = rr_list_gal[i].r;
	  else if(istar < 0.83 * nstar)
	    r_lum = rr_list_gal[i].r;
	}
    }

  for(i = 0; i < num; i++)
    {
      if(P[rr_list_gal[i].index].Type == 0 && r_lum > 0)
	{
	  (*SFRgas) += get_starformation_rate_subfind(rr_list_gal[i].index, &Temperature, &xclouds);
#ifdef DUSTATT
	  if(rr_list_gal[i].r2 > r_lum)
	    rad_bin = DUSTATT - 1;
	  else
	    rad_bin = floor(rr_list_gal[i].r2 / r_lum * (DUSTATT - 2));

	  mcold = xclouds * P[rr_list_gal[i].index].Mass;

	  if(mcold > 0)
	    {
	      Mcold_in_bin[rad_bin] += mcold;
	      Zgas_in_bin[rad_bin] += get_metallicity_subfind(rr_list_gal[i].index) / 0.02 * mcold;
	    }
	}
#endif
    }

#ifdef DUSTATT
  /*- Note that impact parameters > r_lum are not needed, therefore last bin is irrelevant -*/
  for(k = 0; k < DUSTATT - 1; k++)
    {
      if(Mcold_in_bin[k] > 0)
	{
	  Zgas_in_bin[k] /= Mcold_in_bin[k];
	  TauV_in_bin[k] =
	    pow(Zgas_in_bin[k],
		1.6) * (Mcold_in_bin[k] / ((2 * k + 1) * 3.141516 *
					   pow(r_lum / (DUSTATT - 2) * All.Time,
					       2))) * Unit_corrections_for_dust;
	}
      dust[k] = TauV_in_bin[k];
    }
#endif

  for(i = 0, istar = 0; i < num; i++)
    {
      if(P[rr_list_gal[i].index].Type == 4 && r_lum > 0)
	{
	  if(rr_list_gal[i].r <= r_lum)
	    {
	      add_luminosity(rr_list_gal[i].index, &FLAG_AGE, &lum_particle[0], &lum_particle_obs[0]);
	      for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		lum[ilum] += lum_particle[ilum];
#if defined(OBSERVER_FRAME) && !defined(DUSTATT)
	      for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		lum_obs[ilum] += lum_particle_obs[ilum];
#endif
#ifdef DUSTATT
	      rad_bin = floor(rr_list_gal[i].r2 / r_lum * (DUSTATT - 2));
	      rad_bin_z = floor(rr_list_gal[i].r / r_lum * (DUSTATT - 2));

	      if(FLAG_AGE == 1)
		mu = 1;
	      else
		{
		  mu = -1.;
		  while(mu < 0)
		    {
		      mu = gasdev(&mu_seed) * MUWIDTH + MUCENTER;
		      if(mu < 0.1 || mu > 1.0)
			mu = -1.0;
		    }
		}

	      for(k = 0; k < DUSTATT - 1; k++)
		{
		  if(k >= rad_bin_z)
		    term_depth1 += TauV_in_bin[k];
		  if(k <= rad_bin_z)
		    term_depth2 += TauV_in_bin[k];
		  term_depth_total += TauV_in_bin[k];
		}
	      if(term_depth_total > 0.)
		{
		  if(rr_list_gal[i].z >= 0)
		    depth = term_depth1 / (2 * term_depth_total);
		  if(rr_list_gal[i].z < 0)
		    depth = (term_depth2 + term_depth_total) / (2 * term_depth_total);
		}
	      else
		depth = 0.;

	      for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		{
		  Tau_Lambda_over_V = mu * pow(Filters_Effective_Wavelenght[ilum] / 550, -0.7);
		  lum_att[ilum] +=
		    lum_particle[ilum] * exp(-depth * Tau_Lambda_over_V * TauV_in_bin[rad_bin]);
#ifdef OBSERVER_FRAME
		  Tau_Lambda_over_V_obs = mu * pow(Filters_Effective_Wavelenght[ilum] * All.Time / 550, -0.7);
		  lum_obs[ilum] +=
		    lum_particle_obs[ilum] * exp(-depth * Tau_Lambda_over_V_obs * TauV_in_bin[rad_bin]);
#endif
		}
#endif

	      (*meanage) += P[rr_list_gal[i].index].StellarAge;
	      (*meanzstar) += get_metallicity_subfind(rr_list_gal[i].index) / 0.02;
	      istar++;
	    }
	}
    }
  if(istar > 0)
    (*meanage) /= istar;
  if(istar > 0)
    (*meanzstar) /= istar;

  myfree(rr_list_gal);
#endif

#ifdef KD_FRICTION_DYNAMIC
  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      if(P[p].Type == 5)
	{
	  if(*vmax != 0)
	    {
	      BPP(p).BH_sigma = *vmax / sqrt(3);
	      BPP(p).BH_bmax = *halfmassrad;
	    }
	  else
	    {
	      BPP(p).BH_sigma = pow(*totmass / 1e4 / 2.06, 1. / 3.) * 1e3 / sqrt(3);
	      BPP(p).BH_bmax = pow(*totmass * 1e10, 1. / 3.) * 0.023433340 / 2;
	    }
	}
    }
#endif
}


#ifdef LT_ADD_GAL_TO_SUB
void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					       double *pos, double *vel, double *cm, double *veldisp,
					       double *vmax, double *vmaxrad, double *spin,
					       MyIDType * mostboundid, double *halfmassrad, double *mass_tab,
					       double *lum, double *lum_att, double *lum_obs, double *dust,
					       double *meanage, double *meanzstar, double *SFRgas)
#else
void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					       double *pos, double *vel, double *cm, double *veldisp,
					       double *vmax, double *vmaxrad, double *spin,
					       MyIDType * mostboundid, double *halfmassrad, double *mass_tab)
#endif
{
  int i, j, part_index, *npart, numtot, nhalf, offset, num_use, i_use;
  MyIDType mbid;
  double s[3], sloc[3], v[3], vloc[3], max, vel_to_phys, H_of_a, atime;
  double lx, ly, lz, dv[3], dx[3], disp;
  double loclx, locly, loclz, locdisp;
  double boxsize, ddxx;
  sort_r2list *loc_rr_list;
  int minindex, mincpu;
  double mass, mass_tab_loc[6], maxrad, massloc, *masslist;
  MyFloat minpot, *potlist;
  int nstarloc = 0, ndmloc = 0, ngasloc = 0, nbhloc = 0;
  int nstar = 0, ndm = 0, ngas = 0, nbh = 0;
  double rr_tmp, disp_tmp;

#ifdef LT_ADD_GAL_TO_SUB
  float Temperature, xclouds;
  sort_r2list *loc_rr_list_star;
  double lum_particle[LT_ADD_GAL_TO_SUB], lum_particle_obs[LT_ADD_GAL_TO_SUB], r_lum = 0;
  int nlum, ilum, i_star, FLAG_AGE;
  double loc_lum[LT_ADD_GAL_TO_SUB], loc_meanage = 0, loc_meanzstar = 0, loc_SFRgas = 0;
#ifdef OBSERVER_FRAME
  double loc_lum_obs[LT_ADD_GAL_TO_SUB];
#endif

  for(i = 0; i < LT_ADD_GAL_TO_SUB; i++)
    {
      loc_lum[i] = 0;
#ifdef OBSERVER_FRAME
      loc_lum_obs[i] = 0;
#endif
    }
#ifdef DUSTATT
  double loc_lum_att[LT_ADD_GAL_TO_SUB];
  for(i = 0; i < LT_ADD_GAL_TO_SUB; i++)
    loc_lum_att[i] = 0;
  int rad_bin, rad_bin_z, k;
  double mcold, Mcold_in_bin[DUSTATT], Zgas_in_bin[DUSTATT], TauV_in_bin[DUSTATT], total_val;
  double mu, Unit_corrections_for_dust, Tau_Lambda_over_V, rr2_tmp, depth;
#ifdef OBSERVER_FRAME
  double Tau_Lambda_over_V_obs;
#endif
  long mu_seed;
  float MUWIDTH = 0.2, MUCENTER = 0.3;
  float term_depth1 = 0, term_depth2 = 0, term_depth_total = 0;

  for(i = 0; i < DUSTATT; i++)
    {
      Mcold_in_bin[i] = 0;
      Zgas_in_bin[i] = 0;
      TauV_in_bin[i] = 0;
    }

  Unit_corrections_for_dust =
    All.UnitMass_in_g / All.HubbleParam / (1.4 * PROTONMASS * (All.UnitLength_in_cm / All.HubbleParam) *
					   (All.UnitLength_in_cm / All.HubbleParam) * 2.1 * 1.e21);

  mu_seed = -150;

#endif
#endif

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      H_of_a = 0;
      atime = vel_to_phys = 1;
    }

  potlist = (MyFloat *) mymalloc("potlist", NTask * sizeof(MyFloat));

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      if(P[d[i].index].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[d[i].index].u.DM_Potential;
	  minindex = d[i].index;
	}
      switch (P[d[i].index].Type)
	{
	case 0:
	  ngasloc++;
	  break;
	case 1:
	case 2:
	case 3:
	  ndmloc++;
	  break;
	case 4:
	  nstarloc++;
	  break;
	case 5:
	  nbhloc++;
	  break;
	}
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);
  MPI_Allreduce(&ngasloc, &ngas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ndmloc, &ndm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nstarloc, &nstar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nbhloc, &nbh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0, mincpu = -1, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(mincpu < 0)
    {
      printf("ta=%d num=%d\n", ThisTask, num);
      endrun(121);
    }

  if(ThisTask == mincpu)
    {
      for(j = 0; j < 3; j++)
	s[j] = P[minindex].Pos[j];
    }

  MPI_Bcast(&s[0], 3, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);

  /* s[] now holds the position of minimum potential */
  /* we take that as the center */
  for(j = 0; j < 3; j++)
    pos[j] = s[j];


  /* the ID of the most bound particle, we take the minimum binding energy */
  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
#ifdef NO_GAS_CLOUDS
      if(P[d[i].index].Type > 0)
#endif
	if(P[d[i].index].v.DM_BindingEnergy < minpot || minindex == -1)
	  {
	    minpot = P[d[i].index].v.DM_BindingEnergy;
	    minindex = d[i].index;
	  }
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(ThisTask == mincpu)
    {
      if(minindex == -1)	/* This is to cover the quasi impossible case that one cpu have had only gas */
	endrun(875417);		/* particles and the potential in the simulation was larger than 1e30 !      */
      mbid = P[minindex].ID;
    }

  MPI_Bcast(&mbid, sizeof(mbid), MPI_BYTE, mincpu, MPI_COMM_WORLD);

  myfree(potlist);

  *mostboundid = mbid;

  /* let's get bulk velocity and the center-of-mass */

  for(j = 0; j < 3; j++)
    sloc[j] = vloc[j] = 0;

  for(j = 0; j < 6; j++)
    mass_tab_loc[j] = 0;

  for(i = 0, massloc = 0; i < num; i++)
    {
      part_index = d[i].index;

      for(j = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
	  ddxx = P[part_index].Pos[j] - pos[j];
#endif
	  sloc[j] += P[part_index].Mass * ddxx;
	  vloc[j] += P[part_index].Mass * P[part_index].Vel[j];
	}
      massloc += P[part_index].Mass;

      mass_tab_loc[P[part_index].Type] += P[part_index].Mass;
    }

  MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mass_tab_loc, mass_tab, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *totmass = mass;

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];

      s[j] += pos[j];

#ifdef PERIODIC
      while(s[j] < 0)
	s[j] += boxsize;
      while(s[j] >= boxsize)
	s[j] -= boxsize;
#endif

      cm[j] = s[j];
    }


  locdisp = loclx = locly = loclz = 0;


  /* Here we have to perform only on the dm particles for consistency */
  num_use = num;
#ifdef DENSITY_SPLIT_BY_TYPE
  num_use = ndmloc;
#endif

#ifdef LT_ADD_GAL_TO_SUB
  if(nstarloc > 0)
    loc_rr_list_star = (sort_r2list *) mymalloc("loc_rr_list_star", sizeof(sort_r2list) * nstarloc);
  i_star = 0;
#endif

  if(num_use > 0)
    loc_rr_list = (sort_r2list *) mymalloc("loc_rr_list", sizeof(sort_r2list) * num_use);

  for(i = 0, massloc = 0, i_use = 0; i < num; i++)
    {
      part_index = d[i].index;

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
	{
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - s[j]);
#else
	  ddxx = P[part_index].Pos[j] - s[j];
#endif
	  dx[j] = atime * ddxx;
	  dv[j] = vel_to_phys * (P[part_index].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];

	  disp_tmp += P[part_index].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
#ifdef PERIODIC
	  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
	  ddxx = P[part_index].Pos[j] - pos[j];
#endif
	  ddxx = atime * ddxx;
	  rr_tmp += ddxx * ddxx;
	}

      rr_tmp = sqrt(rr_tmp);

#ifdef LT_ADD_GAL_TO_SUB
      if(P[part_index].Type == 4)
	{
	  loc_rr_list_star[i_star].r = rr_tmp;
	  i_star++;
	}
#endif

#ifdef DENSITY_SPLIT_BY_TYPE
      if(P[part_index].Type >= 1 && P[part_index].Type <= 3)  /*-- only for dm part --*/
	{
#endif
	  locdisp += disp_tmp;
	  massloc += P[part_index].Mass;

	  loclx += P[part_index].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
	  locly += P[part_index].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
	  loclz += P[part_index].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

	  loc_rr_list[i_use].mass = P[part_index].Mass;
	  loc_rr_list[i_use].r = rr_tmp;

	  i_use++;
#ifdef DENSITY_SPLIT_BY_TYPE
	}
#endif
    }

  if(i_use != num_use)
    endrun(664321);

#ifdef LT_ADD_GAL_TO_SUB
  if(i_star != nstarloc)
    endrun(774321);
#endif


  MPI_Allreduce(&loclx, &lx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locly, &ly, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loclz, &lz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locdisp, &disp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  npart = (int *) mymalloc("npart", NTask * sizeof(int));

  MPI_Allgather(&num_use, 1, MPI_INT, npart, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, numtot = 0; i < NTask; i++)
    numtot += npart[i];

  if(mass > 0)
    {
      *veldisp = sqrt(disp / (3 * mass));	/* convert to 1d velocity dispersion */

      spin[0] = lx / mass;
      spin[1] = ly / mass;
      spin[2] = lz / mass;

      parallel_sort(loc_rr_list, num_use, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      nhalf = numtot / 2;
      mincpu = 0;

      while(nhalf >= npart[mincpu])
	{
	  nhalf -= npart[mincpu];
	  mincpu++;
	}

/*--- Here we still have to fix for individual masses, 
        maybe we even want the total mass for this ? 
        Note however that even within multi mass DM simulation
        the DM mass within a clean halo should be the same ... ------*/
      if(ThisTask == mincpu)
	*halfmassrad = loc_rr_list[nhalf].r;

      MPI_Bcast(halfmassrad, 1, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);
    }
  else
    *veldisp = *halfmassrad = spin[0] = spin[1] = spin[2] = 0;

  /* compute cumulative mass */

  masslist = (double *) mymalloc("masslist", NTask * sizeof(double));

  for(i = 0, massloc = 0; i < num_use; i++)
    massloc += loc_rr_list[i].mass;

  MPI_Allgather(&massloc, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  for(i = 1; i < NTask; i++)
    masslist[i] += masslist[i - 1];

  for(i = 1; i < num_use; i++)
    loc_rr_list[i].mass += loc_rr_list[i - 1].mass;

  if(ThisTask > 0)
    for(i = 0; i < num_use; i++)
      loc_rr_list[i].mass += masslist[ThisTask - 1];

  for(i = 0, offset = 0; i < ThisTask; i++)
    offset += npart[i];

/*--- Note that here we might want to correct for the baryon fraction ?? ---*/
  for(i = num_use - 1, max = 0, maxrad = 0; i + offset > 5 && i >= 0; i--)
    {
      if(loc_rr_list[i].r <= 0)
	endrun(124523);
      if(loc_rr_list[i].mass / loc_rr_list[i].r > max)
	{
	  max = loc_rr_list[i].mass / loc_rr_list[i].r;
	  maxrad = loc_rr_list[i].r;
	}
    }
  MPI_Allreduce(&max, vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(max < *vmax)
    maxrad = 0;

  *vmax = sqrt(All.G * (*vmax));

  MPI_Allreduce(&maxrad, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  *vmaxrad = max;

#ifdef LT_ADD_GAL_TO_SUB
  if(nstar > 0)
    {

      MPI_Allgather(&nstarloc, 1, MPI_INT, npart, 1, MPI_INT, MPI_COMM_WORLD);

      parallel_sort(loc_rr_list_star, nstarloc, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      nlum = (int) (0.83 * nstar);
      mincpu = 0;

      while(nlum >= npart[mincpu])
	{
	  nlum -= npart[mincpu];
	  mincpu++;
	}

      if(ThisTask == mincpu)
	r_lum = loc_rr_list_star[nlum].r;

      MPI_Bcast(&r_lum, 1, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);

      for(i = 0; i < num; i++)
	{
	  part_index = d[i].index;
	  if(P[part_index].Type == 0 && r_lum > 0)
	    {
	      loc_SFRgas += get_starformation_rate_subfind(part_index, &Temperature, &xclouds);
#ifdef DUSTATT
	      /*- Only 2D projection ! */
	      for(j = 0, rr2_tmp = 0; j < 2; j++)
		{
		  /* for luminosity, take minimum of potential as center */
#ifdef PERIODIC
		  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
		  ddxx = P[part_index].Pos[j] - pos[j];
#endif
		  ddxx = atime * ddxx;
		  rr2_tmp += ddxx * ddxx;
		}

	      rr2_tmp = sqrt(rr2_tmp);

	      if(rr2_tmp > r_lum)
		rad_bin = DUSTATT - 1;
	      else
		rad_bin = floor(rr2_tmp / r_lum * (DUSTATT - 2));

	      mcold = xclouds * P[part_index].Mass;

	      if(mcold > 0)
		{
		  Mcold_in_bin[rad_bin] += mcold;
		  Zgas_in_bin[rad_bin] += get_metallicity_subfind(part_index) / 0.02 * mcold;
		}
#endif
	    }
	}

#ifdef DUSTATT
      /*- Note that impact parameters > r_lum are not needed, therefore last bin is irrelevant -*/
      for(k = 0; k < DUSTATT - 1; k++)
	{
	  MPI_Allreduce(&Mcold_in_bin[k], &total_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  Mcold_in_bin[k] = total_val;

	  if(Mcold_in_bin[k] > 0)
	    {
	      MPI_Allreduce(&Zgas_in_bin[k], &total_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	      Zgas_in_bin[k] = total_val / Mcold_in_bin[k];
	      TauV_in_bin[k] =
		pow(Zgas_in_bin[k],
		    1.6) * (Mcold_in_bin[k] / ((2 * k + 1) * 3.141516 *
					       pow(r_lum / (DUSTATT - 2) * All.Time,
						   2))) * Unit_corrections_for_dust;
	    }

	  dust[k] = TauV_in_bin[k];
	}
#endif

      for(i = 0; i < num; i++)
	{
	  part_index = d[i].index;

	  if(P[part_index].Type == 4 && r_lum > 0)
	    {
	      for(j = 0, rr_tmp = 0; j < 3; j++)
		{
		  /* for luminosity, take minimum of potential as center */
#ifdef PERIODIC
		  ddxx = NEAREST(P[part_index].Pos[j] - pos[j]);
#else
		  ddxx = P[part_index].Pos[j] - pos[j];
#endif
		  ddxx = atime * ddxx;

		  if(j < 2)
		    rr_tmp += ddxx * ddxx;
		}

	      rr_tmp = sqrt(rr_tmp);

	      if(rr_tmp <= r_lum)
		{
		  add_luminosity(part_index, &FLAG_AGE, &lum_particle[0], &lum_particle_obs[0]);
		  for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		    loc_lum[ilum] += lum_particle[ilum];
#if defined(OBSERVER_FRAME) && !defined(DUSTATT)
		  for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		    loc_lum_obs[ilum] += lum_particle_obs[ilum];
#endif

#ifdef DUSTATT
		  rad_bin = floor(rr_tmp / r_lum * (DUSTATT - 2));
		  rad_bin_z = floor(sqrt(ddxx * ddxx + rr_tmp * rr_tmp) / r_lum * (DUSTATT - 2));

		  if(FLAG_AGE == 1)
		    mu = 1;
		  else
		    {
		      mu = -1.;
		      while(mu < 0)
			{
			  mu = gasdev(&mu_seed) * MUWIDTH + MUCENTER;
			  if(mu < 0.1 || mu > 1.0)
			    mu = -1.;
			}
		    }

		  for(k = 0; k < DUSTATT - 1; k++)
		    {
		      if(k >= rad_bin_z)
			term_depth1 += TauV_in_bin[k];
		      if(k <= rad_bin_z)
			term_depth2 += TauV_in_bin[k];
		      term_depth_total += TauV_in_bin[k];
		    }
		  if(term_depth_total > 0.)
		    {
		      if(ddxx >= 0)
			depth = term_depth1 / (2 * term_depth_total);
		      if(ddxx < 0)
			depth = (term_depth2 + term_depth_total) / (2 * term_depth_total);
		    }
		  else
		    depth = 0.;
		  for(ilum = 0; ilum < LT_ADD_GAL_TO_SUB; ilum++)
		    {
		      Tau_Lambda_over_V = mu * pow(Filters_Effective_Wavelenght[ilum] / 550, -0.7);
		      loc_lum_att[ilum] +=
			lum_particle[ilum] * exp(-depth * Tau_Lambda_over_V * TauV_in_bin[rad_bin]);
#ifdef OBSERVER_FRAME
		      Tau_Lambda_over_V_obs =
			mu * pow(Filters_Effective_Wavelenght[ilum] * All.Time / 550, -0.7);
		      loc_lum_obs[ilum] +=
			lum_particle_obs[ilum] * exp(-depth * Tau_Lambda_over_V_obs * TauV_in_bin[rad_bin]);
#endif
		    }
#endif
		}
	      loc_meanage += P[part_index].StellarAge;
	      loc_meanzstar += get_metallicity_subfind(part_index) / 0.02;
	    }
	}
    }

  MPI_Allreduce(loc_lum, lum, LT_ADD_GAL_TO_SUB, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef DUSTATT
  MPI_Allreduce(loc_lum_att, lum_att, LT_ADD_GAL_TO_SUB, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef OBSERVER_FRAME
  MPI_Allreduce(loc_lum_obs, lum_obs, LT_ADD_GAL_TO_SUB, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif

  MPI_Allreduce(&loc_meanage, meanage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_meanzstar, meanzstar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_SFRgas, SFRgas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(nstar > 0)
    {
      (*meanage) /= nstar;
      (*meanzstar) /= nstar;
    }
#endif

#ifdef KD_FRICTION_DYNAMIC
  for(i = 0; i < num; i++)
    {
      part_index = d[i].index;
      if(P[part_index].Type == 5)
	{
	  if(*vmax != 0)
	    {
	      BPP(part_index).BH_sigma = *vmax / sqrt(3);
	      BPP(part_index).BH_bmax = *halfmassrad;
	    }
	  else
	    {
	      BPP(part_index).BH_sigma = pow(*totmass / 1e4 / 2.06, 1. / 3.) * 1e3 / sqrt(3);
	      BPP(part_index).BH_bmax = pow(*totmass * 1e10, 1. / 3.) * 0.023433340 / 2;
	    }
	}
    }
#endif

  myfree(masslist);
  myfree(npart);
  if(num_use > 0)
    myfree(loc_rr_list);

#ifdef LT_ADD_GAL_TO_SUB
  if(nstarloc > 0)
    myfree(loc_rr_list_star);
#endif

}



int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->bound_length > ((struct cand_dat *) b)->bound_length)
    return -1;

  if(((struct cand_dat *) a)->bound_length < ((struct cand_dat *) b)->bound_length)
    return +1;

  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  if(((struct cand_dat *) a)->len > ((struct cand_dat *) b)->len)
    return -1;

  if(((struct cand_dat *) a)->len < ((struct cand_dat *) b)->len)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->subnr < ((struct cand_dat *) b)->subnr)
    return -1;

  if(((struct cand_dat *) a)->subnr > ((struct cand_dat *) b)->subnr)
    return +1;

  return 0;
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran1(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
	*idum = 1;
      else
	*idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ;
	  *idum = IA * (*idum - k * IQ) - IR * k;
	  if(*idum < 0)
	    *idum += IM;
	  if(j < NTAB)
	    iv[j] = *idum;
	}
      iy = iv[0];
    }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}



float gasdev(long *idum)
{
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;

  if(iset == 0)
    {
      do
	{
	  v1 = 2.0 * ran1(idum) - 1.0;
	  v2 = 2.0 * ran1(idum) - 1.0;
	  rsq = v1 * v1 + v2 * v2;
	}
      while(rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0 * log(rsq) / rsq);
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }
  else
    {
      iset = 0;
      return gset;
    }
}

#endif
