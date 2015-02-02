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
 */
package org.biojava.nbio.alignment.aaindex;

import junit.framework.TestCase;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;



public class TestAAINDEXLoading extends TestCase{
/**
 * 
 * M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
 * 
 *   A		R		N		D		C		Q		E		G		H		I		L		K		M		F		P		S		T		W		Y		V
A    2.09
R   -0.50    2.87
N   -0.57    0.60    3.60
D   -0.73    0.13    1.78    4.02
C    0.33   -1.30   -2.08   -2.51    6.99
Q   -0.75    0.13    0.33    0.34   -0.83    2.60
E   -0.12    0.99   -0.16    1.20   -1.97    1.23    2.97
G    0.27   -0.96    0.79   -1.20   -2.11   -0.12   -0.41    4.36
H   -1.42    0.54    0.76   -0.01   -1.50   -0.46   -0.62   -0.40    5.89
I   -0.97   -1.40   -2.43   -2.77    0.13   -1.47   -1.81   -2.93   -1.76    2.76
L   -0.39   -1.19   -2.10   -2.65   -0.31   -1.49   -2.11   -1.98   -0.93    1.56    2.43
K   -0.38    1.42    0.83    0.66   -2.19    0.92    1.11   -0.71    0.31   -1.81   -1.96    2.91
M   -0.04   -0.63   -2.01   -2.58    1.04   -0.13   -1.86   -1.86   -1.04    0.99    1.61   -1.62    3.75
F   -0.76   -1.40   -2.25   -2.19    1.13   -2.31   -1.61   -2.67   -0.22    0.76    1.23   -2.41    0.80    3.28
P   -0.53    0.21   -1.10    0.72   -2.19    0.24   -0.26   -0.04   -1.44   -2.00   -1.56   -0.19   -1.09   -0.91    5.45
S    0.34   -0.06    0.40    0.71    0.31    1.04    0.31    0.29   -0.74   -1.75   -2.30   -0.06   -1.34   -1.11   -0.29    2.36
T    0.13   -0.15    0.30   -0.75   -0.59    0.60   -0.21   -0.81   -0.52   -0.96   -0.86   -0.10   -1.58   -0.69    0.93    1.20    2.04
W   -0.66   -0.04   -2.89   -1.91   -0.76   -0.81   -2.70   -1.21   -1.48    0.25   -0.14   -1.94    0.87    2.29   -5.34   -1.18   -0.57    6.96
Y   -1.25   -0.75   -0.36   -1.21    0.13   -0.61   -1.64   -1.62   -0.12    0.08    0.70   -1.72   -0.41    1.96   -1.98   -1.56   -0.41    2.15    3.95
V    0.02   -1.52   -2.17   -2.02    0.34   -1.38   -1.84   -1.96   -0.35    1.94    0.81   -1.27    0.61    0.51   -1.11   -1.11    0.05   -1.09    0.21    2.05

 */
	public void testSDMmatrix(){
		String matrixName = "PRLA000101";
	
		SubstitutionMatrix<AminoAcidCompound> sdm = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName);
		
		int scale = 1;
		if ( sdm instanceof ScaledSubstitutionMatrix) {
			ScaledSubstitutionMatrix scaledSDM = (ScaledSubstitutionMatrix)sdm;
			scale = scaledSDM.getScale();
			assertEquals(100,scale);
		}
		
		
		AminoAcidCompoundSet aas = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		
		AminoAcidCompound v = aas.getCompoundForString("V");
		AminoAcidCompound w = aas.getCompoundForString("W");
		AminoAcidCompound r = aas.getCompoundForString("R");
		AminoAcidCompound n = aas.getCompoundForString("N");
		
		short rn = sdm.getValue(r,n);		
		assertEquals(60,rn);
		
		short nr = sdm.getValue(n,r);
		assertEquals(rn,nr);
		
		
		short vv = sdm.getValue(v,v); 		
		assertEquals(205,vv);
		
		short vw = sdm.getValue(v,w);
		assertEquals( -109,vw);
		
		short wv = sdm.getValue(w,v);
		assertEquals(vw,wv);
		
	
	
	   
	}
	
	




	/*
	 * 
	 * M rows = ACDEFGHIKLMNPQRSTVWYJ-, cols = ACDEFGHIKLMNPQRSTVWYJ
	 * 
	 * 
	 A		C		D		E		F		G		H		I		K		L		M		N		P		Q		R		S		T		V		W		Y		J 
A   0.224   0.013   0.055   0.068   0.031   0.067   0.048   0.053   0.068   0.050   0.087   0.059   0.067   0.073   0.062   0.074   0.059   0.079   0.033   0.035   0.121
C   0.002   0.739   0.001   0.006   0.012   0.000   0.001   0.004   0.003   0.000   0.000   0.001   0.001   0.005   0.008   0.001   0.001   0.000   0.001   0.000   0.008
D   0.044   0.007   0.284   0.091   0.016   0.041   0.056   0.033   0.034   0.012   0.022   0.094   0.047   0.052   0.025   0.054   0.044   0.025   0.014   0.023   0.030
E   0.052   0.029   0.079   0.251   0.016   0.028   0.026   0.026   0.053   0.019   0.031   0.038   0.037   0.071   0.049   0.031   0.044   0.034   0.010   0.027   0.008
F   0.010   0.029   0.006   0.008   0.291   0.004   0.023   0.046   0.011   0.047   0.032   0.012   0.006   0.010   0.009   0.011   0.013   0.018   0.093   0.073   0.000
G   0.079   0.000   0.066   0.047   0.020   0.455   0.042   0.024   0.033   0.028   0.039   0.073   0.054   0.054   0.040   0.064   0.037   0.039   0.041   0.036   0.038
H   0.013   0.003   0.021   0.011   0.024   0.010   0.284   0.008   0.021   0.011   0.020   0.035   0.008   0.020   0.023   0.013   0.012   0.020   0.014   0.025   0.023
I   0.014   0.016   0.017   0.014   0.058   0.006   0.010   0.235   0.015   0.050   0.048   0.018   0.009   0.009   0.015   0.014   0.023   0.075   0.015   0.030   0.008
K   0.062   0.007   0.039   0.072   0.032   0.027   0.068   0.039   0.294   0.050   0.077   0.055   0.045   0.077   0.122   0.043   0.059   0.044   0.037   0.035   0.053
L   0.028   0.000   0.010   0.017   0.097   0.013   0.024   0.094   0.035   0.311   0.141   0.030   0.030   0.028   0.027   0.019   0.029   0.073   0.064   0.033   0.015
M   0.010   0.000   0.003   0.005   0.015   0.005   0.005   0.020   0.011   0.030   0.167   0.004   0.003   0.017   0.005   0.003   0.007   0.013   0.004   0.008   0.015
N   0.041   0.007   0.080   0.041   0.022   0.044   0.087   0.031   0.042   0.035   0.024   0.239   0.019   0.040   0.031   0.050   0.051   0.021   0.008   0.036   0.030
P   0.053   0.000   0.039   0.036   0.036   0.006   0.027   0.018   0.017   0.034   0.014   0.018   0.412   0.021   0.026   0.037   0.031   0.019   0.018   0.008   0.015
Q   0.040   0.013   0.038   0.060   0.018   0.025   0.042   0.017   0.046   0.026   0.075   0.032   0.019   0.231   0.056   0.032   0.042   0.036   0.007   0.015   0.023
R   0.025   0.023   0.015   0.031   0.010   0.017   0.033   0.017   0.062   0.019   0.015   0.022   0.018   0.047   0.248   0.026   0.028   0.022   0.022   0.023   0.000
S   0.100   0.013   0.088   0.059   0.044   0.073   0.057   0.051   0.062   0.043   0.026   0.096   0.070   0.072   0.079   0.290   0.138   0.057   0.025   0.059   0.053
T   0.054   0.010   0.049   0.058   0.042   0.029   0.037   0.059   0.058   0.039   0.049   0.068   0.042   0.066   0.053   0.099   0.266   0.061   0.021   0.041   0.015
V   0.041   0.000   0.021   0.033   0.040   0.020   0.038   0.148   0.031   0.077   0.051   0.018   0.020   0.043   0.033   0.028   0.044   0.269   0.023   0.049   0.091
W   0.005   0.000   0.002   0.002   0.049   0.006   0.006   0.009   0.006   0.017   0.005   0.003   0.004   0.003   0.008   0.003   0.004   0.003   0.421   0.038   0.000
Y   0.014   0.000   0.013   0.018   0.111   0.012   0.034   0.028   0.017   0.023   0.026   0.023   0.005   0.012   0.024   0.018   0.019   0.028   0.109   0.355   0.023
J   0.002   0.000   0.001   0.000   0.001   0.002   0.002   0.001   0.002   0.001   0.003   0.001   0.001   0.001   0.001   0.002   0.001   0.007   0.000   0.001   0.341
-   0.086   0.092   0.072   0.072   0.045   0.089   0.060   0.043   0.061   0.075   0.048   0.061   0.083   0.050   0.056   0.087   0.050   0.057   0.021   0.048   0.091
//
	 * 
	 * 
	 */
	
	public void testOVEJ920104(){
		String name = "OVEJ920104";
		SubstitutionMatrix<AminoAcidCompound> max = SubstitutionMatrixHelper.getMatrixFromAAINDEX(name);
		
		if ( max instanceof ScaledSubstitutionMatrix) {
			ScaledSubstitutionMatrix scaledMAX = (ScaledSubstitutionMatrix) max;
			int scale = scaledMAX.getScale();
			
			assertEquals(1000,scale);
		}
				
		AminoAcidCompoundSet aas = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		
		AminoAcidCompound minus = aas.getCompoundForString("-");
		
		AminoAcidCompound j = aas.getCompoundForString("J");
		AminoAcidCompound y = aas.getCompoundForString("Y");
		AminoAcidCompound a = aas.getCompoundForString("A");
		
		short ay = max.getValue(a,y);
		assertEquals(35,ay);
		
		
		short ya = max.getValue(y,a);
		assertEquals(14,ya);
		
		short minusa = max.getValue(minus, a);
		assertEquals(86, minusa);
		
		short minusy = max.getValue(minus, y);		
		assertEquals( 48,minusy);
		
		short minusj = max.getValue(minus, j);		
		assertEquals( 91,minusj);
		
	}
	
}
