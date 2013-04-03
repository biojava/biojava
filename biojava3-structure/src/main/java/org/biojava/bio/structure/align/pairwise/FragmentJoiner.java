/*
 *                  BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Jan 28, 2006
 *
 */
package org.biojava.bio.structure.align.pairwise;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.helper.AlignTools;
import org.biojava.bio.structure.align.helper.JointFragments;
import org.biojava.bio.structure.jama.Matrix;


/** Joins the initial Fragments together to larger Fragments
 *
 * @author Andreas Prlic
 * @author Peter Lackner
 * @since 1.5
 * @version %I% %G%
 */
public class FragmentJoiner {

   public static Logger logger =  Logger.getLogger("org.biojava.bio.structure.align");

   public FragmentJoiner() {
      super();

   }


   /**
    * Reallocates an array with a new size, and copies the contents
    * of the old array to the new array.
    * @param oldArray  the old array, to be reallocated.
    * @param newSize   the new array size.
    * @return          A new array with the same contents.
    */
   @SuppressWarnings("rawtypes")
public static Object resizeArray (Object oldArray, int newSize) {
      int oldSize = java.lang.reflect.Array.getLength(oldArray);
      Class elementType = oldArray.getClass().getComponentType();
      Object newArray = java.lang.reflect.Array.newInstance(
            elementType,newSize);
      int preserveLength = Math.min(oldSize,newSize);
      if (preserveLength > 0)
         System.arraycopy (oldArray,0,newArray,0,preserveLength);
      return newArray;
   }


   /**
    *  In helices often many similar fragments can be found. To reduce these to a few
    *  representative ones this check can be used. It does a distance check between
    *  all known Fragments and a new one. If this one is on a similar diagonal and it
    *  has a lower rms, this one is a better representation. Note: shifts of one are
    *  not allowed.
    *
    * @param fragments
    * @param f
    * @param rmsmat
    * @return true - if this is a better representant for a group of locala fragments.
    */
   public static boolean reduceFragments(List<FragmentPair> fragments, FragmentPair f, Matrix rmsmat){
      //logger.info("reducing fragments");
      boolean doNotAdd =false;
      int i = f.getPos1();
      int j = f.getPos2();

      for ( int p =0; p < fragments.size(); p++){
         FragmentPair tmp = (FragmentPair)fragments.get(p);

         int di1 = Math.abs(f.getPos1() - tmp.getPos1());
         int di2 = Math.abs(f.getPos2() - tmp.getPos2());
         //if ((di1 <=4)&& (di2 <=4)){
         if (( Math.abs(di1-di2) == 2)) {
            //System.out.println("mm");
            double rms1 = rmsmat.get(tmp.getPos1(),tmp.getPos2());
            double rms2 = rmsmat.get(i,j);
            doNotAdd = true;
            if ( rms2 < rms1){
               //System.out.println("remplacing " + tmp + " with" + f);
               fragments.remove(p);
               fragments.add(f);
               break;
            }
            p = fragments.size();
         }
      }
      return doNotAdd;
   }


   public JointFragments[] approach_ap3(Atom[] ca1, Atom[] ca2, FragmentPair[] fraglst,
         StrucAligParameters params){

      //System.out  .println("approach ap3");
      //the final list of joined fragments stores as apairs
      List<JointFragments> fll = new ArrayList<JointFragments>();

      FragmentPair[] tmpfidx = new FragmentPair[fraglst.length];
      for ( int i=0 ; i < fraglst.length; i++){
         tmpfidx[i] = (FragmentPair)fraglst[i].clone();
      }

      int n  = tmpfidx.length;

      //boolean[] used = new boolean[n];

      // if two fragments after having been joint have rms < X
      // keep the new fragment.
      //logger.info(" number fragments: " + n);

      for (int i=0;i< fraglst.length;i++){

         /*if ( i%100 == 0 )
                System.out.print(i+ " ");
            if ( i%1000 == 0)
                System.out.println("");
          */
         boolean[] used = new boolean[n];

         int p1i  = tmpfidx[i].getPos1();
         int p1j  = tmpfidx[i].getPos2();
         int l1   = tmpfidx[i].getLength();

         JointFragments f = new JointFragments();

         int maxi=p1i+l1-1;

         f.add(p1i,p1j,0,l1);
         used[i] = true;

         //n = tmpfidx.length;

         for (int j=(i+1);j<n;j++){

            if ( used[j])
               continue;

            int p2i = tmpfidx[j].getPos1();
            int p2j = tmpfidx[j].getPos2();
            int l2  = tmpfidx[j].getLength();

            if ( p2i > maxi) {


               // TODO: replace this with plo angle calculation
               if ( params.isDoAngleCheck()){

                  // was 0.174
                  if (! angleCheckOk(tmpfidx[i], tmpfidx[j], 0.4f))
                     continue;
               }

               if ( params.isDoDistanceCheck()) {

                  if (! distanceCheckOk(tmpfidx[i],tmpfidx[j], params.getFragmentMiniDistance()))
                     continue;
               }

               if ( params.isDoDensityCheck()) {

                  if ( ! densityCheckOk(ca1,ca2, f.getIdxlist(), p2i,p2j, l2, params.getDensityCutoff()))
                     continue;
               }


               if ( params.isDoRMSCheck()) {

                  double rms = rmsCheck(ca1,ca2, f.getIdxlist(), p2i, p2j, l2);
                  if ( rms > params.getJoinRMSCutoff())
                     continue;
                  f.setRms(rms);
               }
               //System.out.println("i: " + i + " j: " + j + " rms: " +  rms);

               //System.out.println("good");

               //System.exit(0);
               f.add(p2i,p2j,0,l2);
               used[j] = true;
               //System.out.println(f);
               maxi = p2i+l2-1;


            }
         }
         fll.add(f);
      }
      //System.out.println("");


      Comparator<JointFragments> comp = new JointFragmentsComparator();
      Collections.sort(fll,comp);
      Collections.reverse(fll);
      int m = Math.min(params.getMaxrefine(),fll.size());
      List<JointFragments> retlst = new ArrayList<JointFragments>();
      for ( int i = 0 ; i < m ; i++){
         JointFragments jf = (JointFragments)fll.get(i);
         //System.out.println(jf);
         retlst.add(jf);
      }

      return (JointFragments[])retlst.toArray(new JointFragments[retlst.size()]);

   }

   private boolean densityCheckOk(Atom[] aa1, Atom[] aa2, List<int[]> idxlist,
         int p2i, int p2j, int l2,
         double densityCutoff){
      JointFragments ftmp = new JointFragments();
      ftmp.setIdxlist(idxlist);
      ftmp.add(p2i,p2j,0,l2);
      AlternativeAlignment ali = new AlternativeAlignment();
      ali.apairs_from_idxlst(ftmp);

      Atom[] aa3 = (Atom[])aa2.clone();

      try {

         int[] idx1 = ali.getIdx1();
         int[] idx2 = ali.getIdx2();

         Atom[] ca1subset = AlignTools.getFragmentFromIdxList(aa1, idx1);

         Atom[] ca2subset = AlignTools.getFragmentFromIdxList(aa3,idx2);

         /*
            ali.calculateSuperpositionByIdx(aa1,aa3);

            Matrix rot = ali.getCurrentRotationMatrix();
            Atom atom = ali.getCurrentShift();

            for ( int i=0;i< ca2subset.length;i++){
                Atom a = ca2subset[i];
                Calc.rotate(a,rot);
                Calc.shift(a,atom);
            }
          */
         double density = getDensity(ca1subset, ca2subset);
         //System.out.println("density " + density);
         if ( density > densityCutoff)
            return false;
         else
            return true;
      } catch (StructureException e){
         e.printStackTrace();
         return false;
      }



   }


   /** this is probably useless
    *
    * @param ca1subset
    * @param ca2subset
    * @return a double
    * @throws StructureException
    */
   private double getDensity(Atom[] ca1subset, Atom[] ca2subset ) throws StructureException{


      //double density = 999;

      Atom centroid1 =  Calc.getCentroid(ca1subset);
      Atom centroid2 = Calc.getCentroid(ca2subset);

      // get Average distance to centroid ...

      double d1 = 0;
      double d2 = 0;

     // double d1n = 0 ;
      //double d2n = 0 ;
      for ( int i = 0 ; i < ca1subset.length;i++){
         double dd1 = Calc.getDistance(centroid1, ca1subset[i]);
         double dd2 = Calc.getDistance(centroid2, ca2subset[i]);

         d1 += dd1;
         d2 += dd2;
         // new
        // for ( int j = 0 ; j < ca1subset.length;j++){
         //   double dd1m = Calc.getDistance(ca1subset[i],ca1subset[j]);
        //    double dd2m = Calc.getDistance(ca2subset[i],ca2subset[j]);
           // d1n += dd1m;
          //  d2n += dd2m;
        // }

      }

      double avd1 = d1 / ca1subset.length;
      double avd2 = d2 / ca2subset.length;

      //double avd1n = Math.sqrt(d1n);
      //double avd2n = Math.sqrt(d2n);


      //if ( Math.min(avd1,avd2)  8)
      //System.out.println(avd1n + " " + avd2n + " avd1 " + avd1 + " avd2 " + avd2);

      //return density;
      //return Math.abs(avd1-avd2);
      return Math.min(avd1,avd2);
   }

   private double rmsCheck(Atom[] a1, Atom[] a2,List<int[]> idxlist, int p2i, int p2j, int l2){

      //System.out.println(dd);
      // check if a joint fragment has low rms ...
      JointFragments ftmp = new JointFragments();
      ftmp.setIdxlist(idxlist);
      ftmp.add(p2i,p2j,0,l2);
      Atom[] a3 = new Atom[a2.length];
      for (int i=0;i < a2.length;i++){
         a3[i] = (Atom)a2[i].clone();
      }
      double rms = getRMS(a1,a3,ftmp);
      return rms;
   }

   /** get the RMS of the JointFragments pair frag
    *
    * @param ca1 the array of all atoms of structure1
    * @param ca2 the array of all atoms of structure1
    * @param frag the JointFragments object that contains the list of identical positions
    * @return the rms
    */
   public static double getRMS(Atom[] ca1, Atom[]ca2,JointFragments frag){
      //      now svd ftmp and check if the rms is < X ...
      AlternativeAlignment ali = new AlternativeAlignment();
      ali.apairs_from_idxlst(frag);
      double rms = 999;
      try {

         int[] idx1 = ali.getIdx1();
         int[] idx2 = ali.getIdx2();

         Atom[] ca1subset = AlignTools.getFragmentFromIdxList(ca1, idx1);

         Atom[] ca2subset = AlignTools.getFragmentFromIdxList(ca2,idx2);

         ali.calculateSuperpositionByIdx(ca1,ca2);

         Matrix rot = ali.getRotationMatrix();
         Atom atom = ali.getShift();

         for ( int i=0;i< ca2subset.length;i++){
            Atom a = ca2subset[i];
            Calc.rotate(a,rot);
            Calc.shift(a,atom);
         }

         rms = SVDSuperimposer.getRMS(ca1subset,ca2subset);
         //double rms = ali.getRMS();
         //System.out.println("rms " + rms);

      } catch (StructureException e ){
         e.printStackTrace();

      }
      return rms;
   }

   public boolean angleCheckOk(FragmentPair a, FragmentPair b, float distcutoff){

      double dist = -999;
      try {
         //Atom v1 = tmpfidx[i].getCenter1();
         //Atom v2 = tmpfidx[j].getCenter2();
         Atom v1 = a.getUnitv();
         Atom v2 = b.getUnitv();
         //System.out.println("v1: "+v1);
         //System.out.println("v2: "+v2);
         dist = Calc.getDistance(v1,v2);

      } catch (StructureException e){
         e.printStackTrace();
         return false;
      }
      if ( dist > distcutoff)
         return false;
      return true;
   }

   private boolean distanceCheckOk(FragmentPair a, FragmentPair b, float fragCompatDist){

      double dd ;
      try {
         Atom c1i = a.getCenter1();
         Atom c1j = b.getCenter1();
         Atom c2i = a.getCenter2();
         Atom c2j = b.getCenter2();
         //System.out.println(c1i + " " + c1j + " " + c2i + " " + c2j);
         dd = Calc.getDistance(c1i,c1j) - Calc.getDistance(c2i,c2j);
         //dd = Calc.getDistance(c1i,c2i) - Calc.getDistance(c1j,c2j);
         //System.out.println(dd);
      } catch (StructureException e){
         e.printStackTrace();
         return false;
      }

      if ( dd < 0)
         dd = -dd;
      if ( dd > fragCompatDist)
         return false;
      return true;

   }

   /**
    * Calculate the pairwise compatibility of fpairs.

     Iterates through a list of fpairs and joins them if
     they have compatible rotation and translation parameters.


    * @param fraglst FragmentPair[] array
    * @param angleDiff angle cutoff
    * @param fragCompatDist distance cutoff
    * @param maxRefine max number of solutions to keep
    * @return JointFragments[]
    */
   public JointFragments[] frag_pairwise_compat(FragmentPair[] fraglst,
         int angleDiff,
         float fragCompatDist,
         int maxRefine  ){


      FragmentPair[] tmpfidx = new FragmentPair[fraglst.length];
      for ( int i=0 ; i < fraglst.length; i++){
         tmpfidx[i] = (FragmentPair)fraglst[i].clone();
      }

      int n  = tmpfidx.length;

      //indicator if a fragment was already used
      int[] used = new int[n];

      //the final list of joined fragments stores as apairs
      List<JointFragments> fll = new ArrayList<JointFragments>();

     // int cnt = 0;
      double adiff = angleDiff * Math.PI / 180d;
      logger.finer("addiff" + adiff);
      //distance between two unit vectors with angle adiff
      double ddiff = Math.sqrt(2.0-2.0*Math.cos(adiff));
      logger.finer("ddiff" + ddiff);

      // the fpairs in the flist have to be sorted with respect to their positions

      while (tmpfidx.length > 0){
         int i = 0;
         int j = 1;
         used[i]=1;
         JointFragments f = new JointFragments();

         //maxi=p1i+tmpfidx[i].l
         //f.add([p1i,p1j],tmpfidx[i].l)

         int p1i  = tmpfidx[i].getPos1();
         int p1j  = tmpfidx[i].getPos2();

         int maxi = p1i+tmpfidx[i].getLength();

         f.add(p1i,p1j,0,tmpfidx[i].getLength());

         n = tmpfidx.length;

         while ( j < n) {
            int p2i = tmpfidx[j].getPos1();
            int p2j = tmpfidx[j].getPos2();
            int l2  = tmpfidx[j].getLength();
            if ( p2i > maxi) {
               try {
                  double dist = Calc.getDistance(tmpfidx[i].getUnitv(), tmpfidx[j].getUnitv());
                  if ( dist  < ddiff) {

                     // centroids have to be closer than fragCompatDist
                     double dd = Calc.getDistance(tmpfidx[i].getCenter1(),tmpfidx[j].getCenter1()) -
                     Calc.getDistance(tmpfidx[i].getCenter2(),tmpfidx[j].getCenter2());
                     if ( dd < 0)
                        dd = -dd;
                     if ( dd < fragCompatDist){
                       // cnt +=1;
                        maxi = p2i+l2;
                        used[j]=1;
                        f.add(p2i,p2j,0,tmpfidx[j].getLength());
                     }
                  }
               } catch (StructureException e) {
                  e.printStackTrace();

               }

            }
            j+=1;
         }

         int red = 0;
         for (int k = 0 ; k < n ; k ++) {
            if (used[k] == 0) {
               tmpfidx[red] = tmpfidx[k];
               red+=1;
            }
         }


         used = new int[n];
         //System.out.println("fragmentjoiner new array size " + red+" old one:"+tmpfidx.length);
         tmpfidx = (FragmentPair[])resizeArray(tmpfidx,red);

         fll.add(f);
      }


      Comparator<JointFragments> comp = new JointFragmentsComparator();
      Collections.sort(fll,comp);
      Collections.reverse(fll);
      int m = Math.min(maxRefine,fll.size());
      List<JointFragments> retlst = new ArrayList<JointFragments>();
      for ( int i = 0 ; i < m ; i++){
         JointFragments jf = (JointFragments)fll.get(i);
         retlst.add(jf);
      }

      return retlst.toArray(new JointFragments[retlst.size()]);

   }


   public void extendFragments(Atom[] ca1, Atom[] ca2 ,JointFragments[] fragments, StrucAligParameters params){

      for(JointFragments p : fragments){
         extendFragments(ca1, ca2, p, params);
      }

   }

   public void extendFragments(Atom[] ca1, Atom[] ca2 , JointFragments fragments, StrucAligParameters params){

      List<int[]> pos = fragments.getIdxlist();

     // int lstart = pos.size();

      // idxlist is always sorted...

      int[] firstP = pos.get(0);
      int pstart1 = firstP[0];
      int pstart2 = firstP[1];

      int[] lastP = pos.get(pos.size()-1);
      int pend1 = lastP[0];
      int pend2 = lastP[1];

      boolean startReached = false;
      boolean endReached   = false;

      while (! (startReached && endReached)){

         if ( ! (startReached) && ((pstart1 <= 0) || (pstart2 <= 0))) {
            startReached = true;

         } else {
            pstart1--;
            pstart2--;
         }

         if ( ! (endReached) && (( pend1 >= (ca1.length -1) ) || ( pend2 >= ca2.length -1  )) ){
            endReached = true;
         } else {
            pend1++;
            pend2++;
         }

         //System.out.println(pstart1 + " " + pstart2 + " " + pend1 + " " + pend2);
         if ( ! startReached){
            double newRms1 = testAdd(ca1, ca2, fragments,pstart1,pstart2);

            if (newRms1 < 3.7 ) {
               fragments.add(pstart1,pstart2,0,1);
            } else {
               startReached = true;
            }
         }

         if( ! endReached){

            double newRms2 = testAdd(ca1, ca2, fragments, pend1, pend2);
            if ( newRms2 < 3.7) {
               fragments.add(pend1,pend2,0,1);
            } else {
               endReached = true;
            }
         }

      }

   }


   private double testAdd(Atom[] ca1, Atom[] ca2, JointFragments fragments, int pstart1, int pstart2){

      JointFragments frag = new JointFragments();
      frag.setIdxlist(fragments.getIdxlist());
      frag.add(pstart1, pstart2, 0,1);

      return FragmentJoiner.getRMS(ca1, ca2, frag);

   }

}




class JointFragmentsComparator
implements Comparator<JointFragments> {

   public int compare(JointFragments one, JointFragments two) {


      int s1 = one.getIdxlist().size();
      int s2 = two.getIdxlist().size();

      double rms1 = one.getRms();
      double rms2 = two.getRms();
      //System.out.print(s1 + " " + s2 + " ");
      if ( s1 > s2 ) {
         //System.out.println("1");
         return 1 ;
      }

      else if ( s1 < s2){
         //System.out.println("-1");
         return -1;
      }
      else{
         if ( rms1 < rms2)
            return 1;
         if ( rms1 > rms2)
            return -1;
         //System.out.println("0");
         return 0;
      }
   }




}




