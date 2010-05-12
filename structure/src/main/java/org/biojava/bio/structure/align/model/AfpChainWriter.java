/*
 *                    BioJava development code
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
 * Created on Feb 15, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.model;

import java.io.StringWriter;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;

public class AfpChainWriter
{

   public static final String newline = System.getProperty("line.separator");

   private static int LINELENGTH = 70;

   public static String toFatCat(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
   {

      boolean printLegend = true;
      boolean longHeader  = true;

      return toFatCatCore(afpChain, ca1, ca2, printLegend, longHeader);   
   }

   public static String toScoresList(AFPChain afpChain){

      // see sippl On distance and similarity in fold space 2008 bioinformatics

      StringWriter writer = new StringWriter();

      writer.append("Sab (nr. equivalent residues): " );
      writer.append(afpChain.getNrEQR()+"");
      writer.append(newline);

      writer.append("Dab (distance between folds a,b): ");
      int dab = afpChain.getCa1Length()+afpChain.getCa2Length() - 2 * afpChain.getNrEQR();
      writer.append(dab+"");
      writer.append(newline);

      writer.append("sab (relative similarity): ");
      double sab = 2 * afpChain.getNrEQR() / (double)( afpChain.getCa1Length() + afpChain.getCa2Length());
      writer.append(sab+"");
      writer.append(newline);

      writer.append("cab (coverage a): ");
      double cab = afpChain.getNrEQR() / (double) afpChain.getCa1Length();
      writer.append(cab+"");
      writer.append(newline);

      writer.append("cba (coverage b): ");
      double cba = afpChain.getNrEQR() / (double) afpChain.getCa2Length();
      writer.append(cba+"");
      writer.append(newline);

      writer.append("seq similarity: ");
      writer.append(afpChain.getSimilarity()+"");
      writer.append(newline);

      writer.append("TM-score: ");
      writer.append(afpChain.getTMScore()+"");
      writer.append(newline);

      return writer.toString();
   }

   public static String toFatCatCore(
         AFPChain afpChain, 
         Atom[] ca1, 
         Atom[] ca2, 
         boolean printLegend, boolean longHeader){

      if(!afpChain.isSequentialAlignment()) {
         //TODO find some topology-independent output format
         return "Can't display circular permutations";
      }

      String name1 = afpChain.getName1();
      String name2 = afpChain.getName2();
      int ca1Length = afpChain.getCa1Length();
      int ca2Length = afpChain.getCa2Length();

      int blockNum = afpChain.getBlockNum();
      int totalLenIni = afpChain.getTotalLenIni();
      double totalRmsdIni = afpChain.getTotalRmsdIni();
      int optLength = afpChain.getOptLength();
      double totalRmsdOpt = afpChain.getTotalRmsdOpt();
      double chainRmsd = afpChain.getChainRmsd();
      double alignScore = afpChain.getAlignScore();
      int alnLength = afpChain.getAlnLength();
      int gapLen = afpChain.getGapLen();
      List<AFP> afpSet = afpChain.getAfpSet();

      double similarity = afpChain.getSimilarity();
      double identity = afpChain.getIdentity();

      if (similarity <0  || identity < 0){
         afpChain.calcSimilarity();
         similarity = afpChain.getSimilarity();
         identity = afpChain.getIdentity();
      }



      String algorithmName = afpChain.getAlgorithmName();
      //String version = afpChain.getVersion();

      double probability = afpChain.getProbability();


      int afpNum = afpSet.size();

      int[] blockGap = afpChain.getBlockGap();


      double[] blockScore = afpChain.getBlockScore();
      double[] blockRmsd = afpChain.getBlockRmsd();
      int[] blockSize = afpChain.getBlockSize();


      int alnbeg1 = afpChain.getAlnbeg1();
      int alnbeg2 = afpChain.getAlnbeg2();

      char[] alnseq1 = afpChain.getAlnseq1();
      char[] alnseq2 = afpChain.getAlnseq2();
      char[] alnsymb = afpChain.getAlnsymb();

      // == end of extractation of data values from afpChain 
      ////////////////////////////////

      StringBuffer txt = new StringBuffer();

      if ( longHeader) {
         txt.append(String.format("Align %s.pdb %d with %s.pdb %d", name1, ca1Length, name2, ca2Length));
      }
      else {
         txt.append(String.format("Align %s.pdb Length1: %d with %s.pdb Length2: %d", name1, ca1Length, name2, ca2Length));
      }
      txt.append(newline);
      if ( afpChain.isShortAlign()){
         txt.append("Short match");
         return txt.toString();
      }
      //txt.append(String.format("raw-score: %.2f norm.-score: %.2f ", alignScore, normAlignScore));

      if ( longHeader ) {
         txt.append(String.format( "Twists %d ini-len %d ini-rmsd %.2f opt-equ %d opt-rmsd %.2f chain-rmsd %.2f Score %.2f align-len %d gaps %d (%.2f%%)",
               blockNum - 1, totalLenIni, totalRmsdIni, optLength, totalRmsdOpt, chainRmsd, alignScore, 
               alnLength, gapLen, (100.0 * (double)gapLen/(double)alnLength)) );
         txt.append(newline);

      }  else {

         if ( ! longHeader)
            printScore(txt,algorithmName,probability,longHeader);

         printScoresInLines(afpChain, blockNum, optLength, totalRmsdOpt, alignScore, alnLength, gapLen, identity, similarity,txt);
      }


      //txt.append(String.format("P-value %.2e Afp-num %d Identity %.2f%% Similarity %.2f%% norm.-score: %.2f"+newline, probability, afpNum, identity * 100, similarity * 100, normAlignScore));

      if ( longHeader) {
         printScore(txt,algorithmName,probability,longHeader);

         txt.append(String.format("Afp-num %d Identity %.2f%% Similarity %.2f%%", afpNum, identity * 100, similarity * 100));
         txt.append(newline);
      } 

      int i;
      double gap;

      if ( longHeader ){
         int fragLen = 8 ; // FatCatParameters.DEFAULT_FRAGLEN;
         for(i = 0; i < blockNum; i ++)  {
            gap = (double)blockGap[i] /( (double)blockGap[i] + fragLen * blockSize[i]);
            txt.append(String.format( "Block %2d afp %2d score %5.2f rmsd %5.2f gap %d (%.2f%%)",
                  i, blockSize[i], blockScore[i], blockRmsd[i], blockGap[i], gap));
            txt.append(newline);
         }
      }

      int     linelen = 70;
      String a;
      String b;
      String c;


      int     t = 0;
      int     ap = alnbeg1;
      int     bp = alnbeg2;
      int     k, len;

      while((alnLength - t) > 0)      {
         if(alnLength - t > linelen)     len = linelen;
         else    len = alnLength - t;


         //System.err.println("t,len:"+t+":"+len);
         a = new String(alnseq1).substring(t,t+len);
         b = new String(alnseq2).substring(t,t+len);
         c = new String(alnsymb).substring(t,t+len);
         //System.err.println("B:" + b);



         txt.append(newline);
         if ( longHeader )
            txt.append(String.format("%14s", " "));
         else 
            txt.append(String.format("%14s", " "));

         if (  longHeader ) {
            for(k = 10; k <= len; k += 10)
               txt.append("    .    :");
            if(k <= len + 5) txt.append("    .");

         } else {

            for(k = 10; k <= len; k += 10)
               txt.append("----+----|");
            if(k <= len + 5) txt.append("----+");


         }

         if ( ap >= ca1.length)
            break;
         if ( bp >- ca2.length)
            break;
         
         String pdb1 = ca1[ap].getParent().getPDBCode();
         String pdb2 = ca2[bp].getParent().getPDBCode();

         txt.append(newline);
         txt.append(String.format("Chain 1:%5s %s"+newline +"%14s%s"+newline+"Chain 2:%5s %s",
               pdb1, a, " ", c, pdb2, b));

         txt.append(newline);
         for(k = 0; k < len; k ++)       {
            if(a.charAt(k) != '-') ap ++;
            if(b.charAt(k) != '-') bp ++;
         }
         t += len;



      }
      txt.append(newline);
      if ( printLegend ){
         if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || 
               algorithmName.equalsIgnoreCase(SmithWaterman3Daligner.algorithmName)){
            txt.append("Note: positions are from PDB; | means alignment of identical amino acids, : of similar amino acids ");

         } else {
            txt.append("Note: positions are from PDB; the numbers between alignments are block index");
         }
         txt.append(newline);
      }
      return txt.toString();

   }

   private static void printScoresInLines(AFPChain afpChain, int blockNum, int optLength, double totalRmsdOpt, double alignScore,
         int alnLength, int gapLen, double identity, double similarity, StringBuffer txt)
   {
      if ( blockNum - 1 > 0) {
         txt.append(String.format( "Twists %d ", blockNum -1 ));
         txt.append(newline);
      }

      txt.append(String.format("Equ: %d ", optLength));
      txt.append(newline);
      txt.append(String.format("RMSD: %.2f ", totalRmsdOpt));
      txt.append(newline);
      txt.append(String.format("Score: %.2f ", alignScore));
      txt.append(newline);
      txt.append(String.format("Align-len: %d ", alnLength));
      txt.append(newline);
      txt.append(String.format("Gaps: %d (%.2f%%)",
            gapLen, (100.0 * (double)gapLen/(double)alnLength)) );
      txt.append(newline);
      if ( afpChain.getTMScore() >= 0) {
         txt.append(String.format("TM-score: %.2f",afpChain.getTMScore()));
         txt.append(newline);
      }
      txt.append(newline);
      txt.append(String.format("Identity: %.2f%% ", identity * 100 ));
      txt.append(newline);
      txt.append(String.format("Similarity: %.2f%%", similarity * 100));
      txt.append(newline);
   }

   private static void printScore(StringBuffer txt,
         String algorithmName,
         double probability, 
         boolean longHeader)
   {
      if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || algorithmName.equalsIgnoreCase(CeSideChainMain.algorithmName) ){
         txt.append(String.format("Z-score %.2f ", probability));
         if ( ! longHeader)
            txt.append(newline);
      } else if ( algorithmName.equalsIgnoreCase(SmithWaterman3Daligner.algorithmName)) {

      } else {
         if ( longHeader ){
            txt.append(String.format("P-value %.2e ",probability));
         }  else {
            txt.append(String.format("P-value: %.2e ",probability));
            txt.append(newline);
         }
      }



   }

   public static String toWebSiteDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2){
      if ( afpChain.getAlgorithmName().equalsIgnoreCase(FatCatFlexible.algorithmName)) {
         String msg =  toFatCat(afpChain,ca1,ca2) ;

         return msg;
      }

      boolean showSeq = true;

      AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2, showSeq);


      //      String msg= toFatCatCore(afpChain,ca1,ca2, printLegend,longHeader);
      //
    

      String msg = toPrettyAlignment(afpChain, ca1, ca2);

            msg = msg + newline + 
            "     | ... Structurally equivalend and identical residues " + newline +
            "     : ... Structurally equivalend and similar residues " + newline + 
            "     . ... Structurally equivalent, but not similar residues. " + newline;

      msg += newline;
      msg += "     To calculate the coordinates of chain 2 aligned on chain 1 apply the following transformation: ";
      msg += newline;
      msg += newline;
      msg += toRotMat(afpChain);
      return msg;


   }

   private static String toPrettyAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
   {
      String name1 = afpChain.getName1();
      String name2 = afpChain.getName2();
      int ca1Length = afpChain.getCa1Length();
      int ca2Length = afpChain.getCa2Length();

      int blockNum = afpChain.getBlockNum();


      int optLength = afpChain.getOptLength();
      double totalRmsdOpt = afpChain.getTotalRmsdOpt();

      double alignScore = afpChain.getAlignScore();
      int alnLength = afpChain.getAlnLength();
      int gapLen = afpChain.getGapLen();


      double similarity = afpChain.getSimilarity();
      double identity = afpChain.getIdentity();

      if (similarity <0  || identity < 0){
         afpChain.calcSimilarity();
         similarity = afpChain.getSimilarity();
         identity = afpChain.getIdentity();
      }


      String algorithmName = afpChain.getAlgorithmName();
      //String version = afpChain.getVersion();

      double probability = afpChain.getProbability();


      // == end of extractation of data values from afpChain 

      StringBuffer txt = new StringBuffer();

      txt.append(String.format("Align %s.pdb Length1: %d with %s.pdb Length2: %d", name1, ca1Length, name2, ca2Length));

      txt.append(newline);

      if ( afpChain.isShortAlign()){
         txt.append("Short match");
         return txt.toString();
      }

      printScore(txt, algorithmName, probability, false);
      printScoresInLines(afpChain, blockNum, optLength, totalRmsdOpt, alignScore, alnLength, gapLen,identity, similarity, txt);
      txt.append(newline);

      int[] optLen = afpChain.getOptLen();
      int[][][] optAln = afpChain.getOptAln();


      int i, j,p1, p2;

      int k;
      int p1b = 0;
      int p2b = 0;

      int     len = 0;
      StringWriter alnseq1 = new StringWriter();
      StringWriter alnseq2 = new StringWriter();
      StringWriter alnsymb = new StringWriter();
      StringWriter header1  = new StringWriter();
      StringWriter footer1  = new StringWriter();
      StringWriter header2  = new StringWriter();
      StringWriter footer2  = new StringWriter();
      StringWriter block    = new StringWriter();

      for(i = 0; i < blockNum; i ++)  {   

         for(j = 0; j < optLen[i]; j ++) {

            p1 = optAln[i][0][j];
            p2 = optAln[i][1][j];

//               System.out.println(p1 + " " + p2 + " " +  footer2.toString());

            if ( len == 0){
               //the first position of sequence in alignment
               formatStartingText(p1,p2,header1,header2,footer1,footer2,ca1,ca2);
            } else {
               // check for gapped region
               int lmax = (p1 - p1b - 1)>(p2 - p2b - 1)?(p1 - p1b - 1):(p2 - p2b - 1);
               for(k = 0; k < lmax; k ++)      {

                  
                  formatGappedRegion(ca1, ca2, txt, p1, p2, k, p1b, p2b, alnseq1, alnseq2, alnsymb, header1, footer1, header2,
                        footer2, block,len);           
                  len++;
                  doLenCheck(len,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block)  ;              
               }
            }
          
            // ALIGNED REGION
//           System.out.println(len + " >" + header1.toString() +"< ");
//           System.out.println(len + " >" + header2.toString() +"< ");   
//           System.out.println(len + " >" + alnseq1.toString() +"< ");
//           System.out.println(len + " >" + alnsymb.toString() +"< ");
//           System.out.println(len + " >" + alnseq2.toString() +"< ");
//           System.out.println(len + " >" + footer1.toString() +"< ");
            formatAlignedRegion(ca1, ca2, p1, p2, alnseq1, alnseq2, alnsymb, header1, footer1, header2, footer2, block,len);
//            System.out.println(len + " >" + header1.toString() +"< ");
//            System.out.println(len + " >" + header2.toString() +"< ");   
//            System.out.println(len + " >" + alnseq1.toString() +"< "); 
//            System.out.println(len + " >" + alnsymb.toString() +"< ");
//            System.out.println(len + " >" + alnseq2.toString() +"< ");
//            System.out.println(len + " >" + footer1.toString() +"< ");
            
            len++;
            
            doLenCheck(len,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block)  ;

            p1b = p1;
            p2b = p2;

            //header1.append(newline);
            //header2.append(newline);

         }

      }

      alnLength = len;

      doLenCheck(LINELENGTH,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block);
      return txt.toString();
   }

   private static void formatGappedRegion(Atom[] ca1, Atom[] ca2, StringBuffer txt, int p1, int p2, int k, int p1b, int p2b,
         StringWriter alnseq1, StringWriter alnseq2, StringWriter alnsymb, StringWriter header1, StringWriter footer1,
         StringWriter header2, StringWriter footer2, StringWriter block, int len)
   {

      // DEAL WITH GAPS
      int tmppos = (p1 - p1b - 1);
      block.append("g");

      if(k >= tmppos) {
         //alnseq1[len] = '-';
         alnseq1.append('-'); 
         header1.append(" ");
         header2.append(" ");

      }
      else {
         int  pos1=p1b+1+k ;
         char oneletter = getOneLetter(ca1[pos1].getParent());
         alnseq1.append(oneletter);

         formatPosition(pos1,ca1, len, header1, header2);
                
      }

      if(k >= (p2 - p2b - 1)) {
         //alnseq2[len] = '-';
         alnseq2.append('-');
         footer1.append(" ");
         footer2.append(" ");

      }
      else  {
         int pos2=p2b+1+k;
         char oneletter = getOneLetter(ca2[pos2].getParent());

         alnseq2.append(oneletter);

         formatPosition(pos2, ca2, len, footer1, footer2);
        
      }
      //alnsymb[len ++] = ' ';
      alnsymb.append(' ');

   }

   private static Integer getPDBResnum(Group g){
      Integer pos1 = null;
      if ( g != null ){
         String pdbResNum = g.getPDBCode();
         try {
        
            pos1 = Integer.parseInt(pdbResNum);
           
         } catch (NumberFormatException e){
            // can be ignored. this can be the case with insertion codes.           
         }
      }
      return pos1;
   }
   
   private static void formatPosition(int pos1, Atom[] ca, int len, StringWriter header1, StringWriter header2)
   {
      int linePos = len % LINELENGTH;
  
      if ( header1.getBuffer().length() < linePos) {
         // fill up the buffer, we are probably shortly after the start...
         for ( int i = header1.getBuffer().length() ; i< linePos ; i++){
            header1.append(" ");
         }
      }
  
      
      
      Atom a = ca[pos1];
      Group g = a.getParent();
      
      boolean hasInsertionCode = false;
      
      Integer pdbPos = getPDBResnum(g);
      if ( pdbPos == null) {
         hasInsertionCode = true;
      } else {
         pos1= pdbPos.intValue();
      }
      
      if ( (pos1 %10  == 0) && ( ! hasInsertionCode)) {
         CharSequence display = getPDBPos(a);
         
         boolean ignoreH1 = false; 
         
         // make sure we don't have a problem with the left boundary...
         if ( header1.getBuffer().length()-1 > linePos) {
            ignoreH1 = true;
            System.out.println("Ignore h1: " + len + " " + header1.getBuffer().length() + " linePos: " + linePos +"  >" + header1.toString() +"<");
         }
         //System.out.println(len + " p1:" + tmp + " = " + pos1 + " " + " " + display + " " + ignoreH1);
         if ( ! ignoreH1) {
            header1.append(String.format("%-10s",display ));
            header2.append("|");
         } else {
            header2.append("|");
         }
        
      } else if ( hasInsertionCode){
         String insCode = getInsertionCode(g);
         if ( insCode.length() == 1)
            header2.append(insCode);
         else {
            header2.append("!");
         }
      } else if ( ((pos1) %5 ) == 0 && len > 5) {
         header2.append(".");
      } else {
         if ( len > 0)
            header2.append(" ");
      }
      
   }

   private static String getInsertionCode(Group g)
   {
      
      String pdbResNum = g.getPDBCode();
      StringBuffer strBuff = new StringBuffer();
      char c;
      
      for (int i = 0; i < pdbResNum.length() ; i++) {
          c = pdbResNum.charAt(i);
          
          if (! Character.isDigit(c)) {
              strBuff.append(c);
          }
      }
      return strBuff.toString();
      
   }

   private static void formatAlignedRegion(Atom[] ca1, Atom[] ca2, int p1, int p2, StringWriter alnseq1, StringWriter alnseq2,
         StringWriter alnsymb, StringWriter header1, StringWriter footer1, StringWriter header2, StringWriter footer2, StringWriter block, int len)
   {
      char c1 =  getOneLetter(ca1[p1].getParent());
      char c2 =  getOneLetter(ca2[p2].getParent());

      alnseq1.append(c1);              
      alnseq2.append(c2);


      if ( c1 == c2){               
         alnsymb.append('|');
         //alnsymb[len ++] = '|';
      } else {
         double score = AFPAlignmentDisplay.aaScore(c1,c2);

         if ( score > 1)
            alnsymb.append( ':');
         else 
            alnsymb.append( '.');               
      }


      formatPosition(p1, ca1,len, header1, header2);

      formatPosition(p2,ca2,len, footer1, footer2);

   }

   private static void formatStartingText(int p1, int p2, StringWriter header1, StringWriter header2, StringWriter footer1,
         StringWriter footer2, Atom[] ca1, Atom[] ca2)
   {

      header1.append(String.format("%-10s", getPDBPos(ca1[p1])));
      header2.append("|");
      footer1.append(String.format("%-10s", getPDBPos(ca2[p2])));
      footer2.append("|");
//
//     
//
//      Integer pdb1 = getPDBResnum(ca1[p1].getParent());
//      Integer pdb2 = getPDBResnum(ca2[p2].getParent());
//      
//      if ( pdb1 != null)
//         p1 = pdb1.intValue();
//      if ( pdb2 != null)
//         p2 = pdb2.intValue();
      
//      int left1 = p1 % 10;
//      if ( left1 < 5) {
//         int space1 = 11 - left1 ;
//         String f1 = "%-" + space1 + "s";
//         header1.append(String.format(f1,getPDBPos(ca1[p1])));
//         //System.out.println(">"+header1.toString()+"<" + f1);
//      } else {
//         for ( ; left1 < 10 ; left1++){
//            header1.append("-");
//         }
//      }
      
      
//      int left2 = p2 % 10;
//      if ( left2 < 5 ) {
//         int space2 = 11 - left2 ;
//         String f2 = "%-" + space2 + "s";
//         footer1.append(String.format(f2,getPDBPos(ca2[p2])));
//      } else {
//         for ( ; left2 < 10 ; left2++){
//            footer1.append("-");
//         }
//      }


      //System.out.println("start at:" + p1 + " " + left1 + " " + p2 + " " + left2);
      //      if ( p1 > 0 ){
      //         header1.append(String.format("%-10s",getPDBPos(ca1[p1])));
      //         header2.append("|");
      //      }
      //      if ( p2 > 0 ) {
      //         footer1.append(String.format("%-10s",getPDBPos(ca2[p2])));
      //
      //         footer2.append("|");
      //      }

   }

   private static boolean doLenCheck(int len, StringBuffer txt, StringWriter header1, StringWriter header2, StringWriter alnseq1,
         StringWriter alnsymb, StringWriter alnseq2, StringWriter footer1, StringWriter footer2, StringWriter block)
   {

      if ( len % LINELENGTH  == 0) {

         //txt.append("|");
         txt.append(header1);
         //txt.append("|");
         txt.append(newline);
         //txt.append("|");
         txt.append(header2);
         //txt.append("|");
         txt.append(newline);
         //txt.append("|");
         txt.append(alnseq1);
         //txt.append("|");
         txt.append(newline);

         //txt.append("|");
         txt.append(alnsymb);
//         txt.append(newline);
//         txt.append(block);
         //txt.append("|");
         txt.append(newline);
         //txt.append("|");
         txt.append(alnseq2);
         //txt.append("|");
         txt.append(newline);
         //txt.append("|");
         txt.append(footer2);
         //txt.append("|");
         txt.append(newline);
         //txt.append("|");
         txt.append(footer1);
         //txt.append("|");
         txt.append(newline);
         txt.append(newline);
         txt.append(newline);


         alnseq1.getBuffer().replace(0, LINELENGTH, "");
         alnseq2.getBuffer().replace(0, LINELENGTH, "");         
         alnsymb.getBuffer().replace(0, LINELENGTH, "");              
         header1.getBuffer().replace(0, LINELENGTH, "");
         header2.getBuffer().replace(0, LINELENGTH , "");
         footer1.getBuffer().replace(0, LINELENGTH, "");         
         footer2.getBuffer().replace(0, LINELENGTH, "");
         block.getBuffer().replace(0, LINELENGTH, "");

         StringBuffer buf = header1.getBuffer();
         for ( int i=0;i<buf.length();i++){
            char c = buf.charAt(i);
            if ( c != ' '){
               buf.setCharAt(i, ' ');
            }
         }
         buf = footer1.getBuffer();
         for ( int i=0;i<buf.length();i++){
            char c = buf.charAt(i);
            if ( c != ' '){
               buf.setCharAt(i, ' ');
            }
         }

         return true;
      }

      return false;


   }

   private static CharSequence getPDBPos(Atom atom)
   {

      Group g = atom.getParent();
      if ( g!= null){
         Chain c = g.getParent();
         if (c != null){
            return g.getPDBCode()+":" + c.getName() ;
            //return g.getPDBCode()+":" + c.getName() + "." + getOneLetter(g) ; 
         }
      }
      return "!";
   }

   private static char getOneLetter(Group g){

      try {
         Character c = StructureTools.get1LetterCode(g.getPDBName());
         return c;
      } catch (Exception e){
         return 'X';
      }
   }


   public static String toDBSearchResult(AFPChain afpChain)
   {
      StringBuffer str = new StringBuffer();

      str.append(afpChain.getName1());
      str.append("\t");
      str.append(afpChain.getName2());
      str.append("\t");
      str.append(String.format("%.2f",afpChain.getAlignScore()));
      str.append("\t");     
      if ( afpChain.getAlgorithmName().equalsIgnoreCase(CeMain.algorithmName)){
         str.append(String.format("%.2f",afpChain.getProbability()));
      } else {
         str.append(String.format("%.2e",afpChain.getProbability()));
      }
      str.append("\t");
      str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
      str.append("\t");
      str.append(afpChain.getCa1Length());
      str.append("\t");
      str.append(afpChain.getCa2Length());      
      str.append("\t");
      str.append(afpChain.getSimilarity1());
      str.append("\t");
      str.append(afpChain.getSimilarity2());
      str.append("\t");
      str.append(newline);

      return str.toString();
   }

   public static String toRotMat(AFPChain afpChain)
   {

      Matrix[] blockRotationMatrix = afpChain.getBlockRotationMatrix();
      int blockNum = afpChain.getBlockNum();
      Atom[] blockShiftVector = afpChain.getBlockShiftVector();

      StringBuffer txt = new StringBuffer();

      if ( blockRotationMatrix == null || blockRotationMatrix.length < 1)
         return "";


      for ( int blockNr = 0 ; blockNr < blockNum  ; blockNr++){
         Matrix m = blockRotationMatrix[blockNr];
         Atom shift   = blockShiftVector[blockNr];
         if ( blockNum > 1) {
            txt.append("Operations for block " );
            txt.append(blockNr);
            txt.append(newline);
         }

         String origString = "orig";
         if ( blockNr > 0)
            origString = (blockNr)+""; 


         txt.append(String.format("     X"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,0),m.get(1,0), m.get(2,0), shift.getX()));
         txt.append( newline); 
         txt.append(String.format("     Y"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,1),m.get(1,1), m.get(2,1), shift.getY()));
         txt.append( newline);
         txt.append(String.format("     Z"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,2),m.get(1,2), m.get(2,2), shift.getZ()));
         txt.append(newline);
      }
      return txt.toString();
   }

   public static String toCE(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
   {



      String name1 = afpChain.getName1();
      String name2 = afpChain.getName2();

      int optLength = afpChain.getOptLength();
      double totalRmsdOpt = afpChain.getTotalRmsdOpt();

      int alnLength = afpChain.getAlnLength();
      int gapLen = afpChain.getGapLen();


      double similarity = afpChain.getSimilarity();
      double identity = afpChain.getIdentity();
      if (similarity == -1 || identity == -1){
         afpChain.calcSimilarity();
         similarity = afpChain.getSimilarity();
         identity = afpChain.getIdentity();
      }


      double probability = afpChain.getProbability();


      int alnbeg1 = afpChain.getAlnbeg1();
      int alnbeg2 = afpChain.getAlnbeg2();

      char[] alnseq1 = afpChain.getAlnseq1();
      char[] alnseq2 = afpChain.getAlnseq2();


      long calculationTime = afpChain.getCalculationTime();

      // == end of extractation of data values from afpChain 



      StringBuffer txt = new StringBuffer();

      txt.append("Chain 1: ");
      txt.append(name1);
      txt.append(" (Size=");
      txt.append(ca1.length);
      txt.append(")");
      txt.append(newline);
      txt.append("Chain 2: ");
      txt.append(name2);
      txt.append(" (Size=");
      txt.append(ca2.length);
      txt.append(")");
      txt.append(newline);
      txt.append(newline);
      txt.append(String.format("Alignment length = %d Rmsd = %.2fA Z-Score = %.1f",optLength,totalRmsdOpt,probability));
      txt.append(String.format(" Gaps = %d(%.1f%%) CPU = %d ms. Sequence identities = %.1f%%",gapLen,( gapLen*100.0/optLength),calculationTime,identity*100));

      int     linelen = 70;
      String a;
      String b;



      int     t = 0;
      int     ap = alnbeg1;
      int     bp = alnbeg2;
      int     k, len;

      while((alnLength - t) > 0)      {
         if(alnLength - t > linelen)     len = linelen;
         else    len = alnLength - t;


         //System.err.println("t,len:"+t+":"+len);
         a = new String(alnseq1).substring(t,t+len);
         b = new String(alnseq2).substring(t,t+len);

         //System.err.println("B:" + b);

         /*
            txt.append(newline);
            txt.append(String.format("%14s", " "));

            for(k = 10; k <= len; k += 10)
                txt.append("    .    :");
            if(k <= len + 5) txt.append("    .");
          */

         //String pdb1 = ca1[ap].getParent().getPDBCode();
         //String pdb2 = ca2[bp].getParent().getPDBCode();
         txt.append(newline);
         txt.append(String.format("Chain 1:%5s %s"+newline+"Chain 2:%5s %s",
               (ap+1), a, (bp+1), b));
         txt.append(newline);
         for(k = 0; k < len; k ++)       {
            if(a.charAt(k) != '-') ap ++;
            if(b.charAt(k) != '-') bp ++;
         }
         t += len;

      }
      txt.append(newline);

      txt.append(toRotMat(afpChain));

      return txt.toString();


   }



}
