/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
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
 *
 * Created on Jun 17, 2009
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */

package org.biojava.bio.structure.align.fatcat.calc;

import org.biojava.bio.structure.align.model.AFPChain;


/* CDF: cumultive density function, SF: survial function for EVD distribution, Gumbel maximum function*/

/* EVD (extreme value distribution)                                                                                                                                                                                    
 * Gumbel maximum distribution (longer right tail) (Type II)                                                                                                                                                           
 * pdf (probability density function) f(x) = 1/b exp(-(x-u)/b)) * exp(-exp(-(x-u)/b)))                                                                                                                                 
 * cdf (cumulative distribution function)  F(x) = exp(-exp(-x))                                                                                                                                                        
 * sf (survial function)              S(x) = 1 - F(x) = 1 - exp(-exp(-x))                                                                                                                                              
 * Gumbel minimum distribution (longer left tail) (Type I)                                                                                                                                                             
 * pdf (probability density function) f(x) = 1/b exp((x-u)/b)) * exp(-exp((x-u)/b)))                                                                                                                                   
 * cdf (cumulative distribution function)  F(x) = 1 - exp(-exp(x))                                                                                                                                                     
 * sf (survial function)              S(x) = 1 - F(x) = exp(-exp(x))                                                                                                                                                   
 * the FATCAT follows Gumbel maximum distribution                                                                                                                                                                      
 * The survival function is the probability that the variate takes a value greater than x.                                                                                                                             
 */
public class SigEva
{


   double  mu;
   double  beta;
   double  mu_a;
   double  beta_a;
   double  mu_b;
   double  beta_b;

   public SigEva(){
      getPara(3,0);
   }


   public void getPara(int set, int len)
   {
      if(set == 1)    {
         mu_a = 0.2461;
         mu_b = 17.1530;
         beta_a = 0.1284;
         beta_b = 1.3756;
      } //fitting, based on the old benchmark (the fold.list is redundant), old parameters                                                                                                                           
      else if(set == 2)       {
         mu_a = 1.1137;
         mu_b = 6.5574;
         beta_a = 0.6448;
         beta_b = -12.2793;
      } //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, "NS3"                                                                                                                               
      else if(set == 3)       {
         mu_a = 0.8440;
         mu_b = 30.2160;
         beta_a = 0.3525;
         beta_b = 18.6652;
      } //score * sqrt(optLen / (rmsd * (twist + 1))), the new best parameters for flexible FATCAT, "FNS8"                                                                                                           
      else if(set == 4)       {
         mu_a = 0.4708;
         mu_b = 38.9863;
         beta_a = 0.2511;
         beta_b = 12.9228;
      } //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=2, modified scoring                                                                                                 
      else if(set == 5)       {
         mu_a = 1.3794;
         mu_b = -14.4778;
         beta_a = 0.7465;
         beta_b = -22.9452;
      } //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=2, modified scoring 3                                                                                            
      else if(set == 6)       {
         mu_a = 0.6036;
         mu_b = 35.4783;
         beta_a = 0.3136;
         beta_b = 13.3922;
      } //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=1                                                                                                                   
      else if(set == 7)       {
         mu_a = 0.7183;
         mu_b = 27.9647;
         beta_a = 0.4688;
         beta_b = -1.3293;
      } //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=1                                                                                                                
      else if(set == 8)       {
         mu_a = 0.4813;
         mu_b = 34.2051;
         beta_a = 0.2618;
         beta_b = 12.4581;
      } //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=3, modified scoring                                                                                                 
      else if(set == 9)       {
         mu_a = 0.6672;
         mu_b = 26.5767;
         beta_a = 0.4373;
         beta_b = -1.4017;
      } //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=3, modified scoring 2                                                                                            
      else    {
         System.err.println("no corresponding parameter set found!");
      }

      mu = beta = .0;
      calMu(len);
      calBeta(len);
   }


   private void calMu(int len)
   {
      mu = mu_a * len + mu_b;
      //printf("mu = %.4f\n", mu);                                                                                                                                                                                   
   }

   private void calBeta(int len)
   {
      beta = beta_a * len + beta_b;
      //printf("beta = %.4f\n", beta);                                                                                                                                                                               
   }


   private int aveLen(int len1, int len2)
   {
      //int   len = (int(sqrt(len1 * len2)));                                                                                                                                                                        
      //int   len = len1 < len2?len1:len2; //use the minimum length, bad discriment                                                                                                                                  
      int     len = (int)(0.5 * (len1 + len2));
      return len;
   }

   public double calSigAll(FatCatParameters params, AFPChain afpChain){
      
      int twist = params.getMaxTra();
      int sparse = params.getSparse();
      int len1 = afpChain.getCa1Length();
      int len2 = afpChain.getCa2Length();
      double score = afpChain.getAlignScore();
      double rmsd = afpChain.getTotalRmsdOpt();
      int optLen = afpChain.getOptLength();
      int r = afpChain.getBlockNum() -1;
      
      return calSigAll(twist,sparse,len1, len2,score,rmsd,optLen,r);
      
   }
   
   private double calSigAll(int twist, int sparse, int len1, int len2, double score, double rmsd, int optLen, int r)
   {
      int     len = aveLen(len1, len2);
      if(sparse == 2) { //use sparse sampling = 2                                                                                                                                                                    
         if(twist == 0)  getPara(4, len); //rigid-FATCAT                                                                                                                                                        
         else    getPara(5, len);        //flexible-FATCAT                                                                                                                                                      
      }
      else if(sparse == 3)    { //sparse sampling = 3                                                                                                                                                                
         if(twist == 0)  getPara(8, len); //rigid-FATCAT                                                                                                                                                        
         else    getPara(8, len); //flexible-FATCAT                                                                                                                                                             
      }
      else if(sparse == 1)    { //sparse sampling = 1                                                                                                                                                                
         if(twist == 0)  getPara(6, len); //rigid-FATCAT                                                                                                                                                        
         else    getPara(7, len); //flexible-FATCAT                                                                                                                                                             
      }
      else    { //no sparse sampling or the corresponding sparse parameters are not fitted                                                                                                                           
         if(twist == 0)  getPara(2, len); //rigid-FATCAT                                                                                                                                                        
         else    getPara(3, len);        //flexible-FATCAT                                                                                                                                                      
      }

      double  mods = normScore(score, rmsd, optLen, r);

      double  t = (mods - mu) / beta;

      double  sf = sF(t);
      return sf;
   }

   private double sF(double t){
      return (1- cDF(t));
   }
   
   private double cDF(double t){
      return (Math.exp(-Math.exp(-t)));
   }

   //the chaining score is normalized by rmsd, twist and optimal alignment length                                                                                                                                         
   private double  normScore(double score, double rmsd, int optLen, int r)
   {
      //double        score1 = modScore(score, r);                                                                                                                                                                   
      double  score1 = score;
      if(r > 0)       score1 /= Math.sqrt(r + 1);
      //it is tested that flexible score is more linear relevant to 1/r2 than 1/r                                                                                                                                   
      if(rmsd < 0.5)  score1 *= Math.sqrt((optLen) / 0.5);
      else    score1 *= Math.sqrt((optLen) / rmsd);
      return score1;
   }

   
   public double calNS(FatCatParameters params, AFPChain afpChain){
      int len1 = afpChain.getCa1Length();
      int len2 = afpChain.getCa2Length();
      double score =afpChain.getAlignScore();
      double rmsd = afpChain.getTotalRmsdOpt();
      int optLen  = afpChain.getOptLength();
      int r = afpChain.getBlockNum() -1;
      return calNS(len1, len2,score,rmsd,optLen,r);
   }
   private double calNS(int len1, int len2, double score, double rmsd, int optlen, int r)
   {
      int     len = aveLen(len1, len2);
      getPara(3, len);
      double  ns0 = normScore(score, rmsd, optlen, r);
      double  ns1 = (ns0 - mu) / beta;
      return ns1;
   }
}
