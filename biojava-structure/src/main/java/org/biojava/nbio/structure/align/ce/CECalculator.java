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
 * Created on Sep 25, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.alignment.aaindex.ScaledSubstitutionMatrix;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFP;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;

import java.util.ArrayList;
import java.util.List;



/** This is based on the original Combinatorial Extension (CE) source code from 2003 or 2004 (CE version 2.3),
 * as has been originally developed by I. Shindyalov and P.Bourne (1998).
 * The original CE paper is available from here: <a href="http://peds.oxfordjournals.org/cgi/content/short/11/9/739">http://peds.oxfordjournals.org/cgi/content/short/11/9/739</a>.
 * 
 * This class is a pretty much exact 1:1 port from C, where I cared about exact reproduce of the CE results
 * and not about Java style.
 * 
 * @author Andreas Prlic

 *
 */
public class CECalculator {

	protected static final boolean isPrint = false;
	private static final boolean debug = false;  

	int[] f1;
	int[] f2;
	double[][]dist1;
	double[][]dist2;
	protected double[][]mat;
	protected int[] bestTrace1;
	protected int[] bestTrace2;
	protected int[][] bestTraces1;
	protected int[][] bestTraces2;
	protected int nBestTrace;
	protected int nBestTraces;
	double d_[] = new double[20];
	protected int[] bestTracesN;
	protected double bestTraceScore;
	protected int nTrace;
	protected double[] bestTracesScores;
	protected int[] trace1;
	protected int[] trace2;

	protected static final 	double  zThr=-0.1;

	long timeStart ;
	long timeEnd;
	private int nAtom;

	// the equivalent positions in the alignment...
	private int[] align_se1;
	private int[] align_se2;


	private int lcmp;
	private int[] bestTraceLen;
	private Matrix r;
	private Atom t;
	protected int nTraces;

	private double z;
	private int[] traceIndexContainer;

	protected CeParameters params;
	// SHOULD these fields be PARAMETERS?

	protected static final int nIter = 1;
	private static final boolean distAll = false;

	List<MatrixListener> matrixListeners;



	public CECalculator(CeParameters params){
		timeStart = System.currentTimeMillis();
		dist1= new double[0][0];
		dist2= new double[0][0];
		this.params = params;
		matrixListeners = new ArrayList<MatrixListener>();

	}

	/**
	 * 
	 * @param afpChain A new AFPChain, which will be filled in by this function
	 * @param ca1
	 * @param ca2
	 * @return afpChain
	 * @throws StructureException
	 */
	public AFPChain extractFragments(AFPChain afpChain,
			Atom[] ca1, Atom[] ca2) throws StructureException{

		int nse1 = ca1.length;
		int nse2 = ca2.length;

		afpChain.setCa1Length(nse1);
		afpChain.setCa2Length(nse2);

		int traceMaxSize=nse1<nse2?nse1:nse2;

		f1 = new int[nse1];
		f2 = new int[nse2];
		
		dist1 = initIntraDistmatrix(ca1, nse1);
		dist2 = initIntraDistmatrix(ca2, nse2);
		
		
		if ( debug )
			System.out.println("parameters: " + params);
		
		if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.SEQUENCE_CONSERVATION){
			if ( params.getSeqWeight() < 1)
				params.setSeqWeight(2);
		}
		
		int winSize = params.getWinSize();

		int winSizeComb1 = (winSize-1)*(winSize-2)/2;		

		traceIndexContainer = new int[traceMaxSize];

		// CE: unused code. distAll is always false and both loops do the same???
		// CE v2.3 calls this Weight factors for trace extension
		if(distAll ) {
			for(int i=0; i<traceMaxSize; i++)
				traceIndexContainer[i]=(i+1)*i*winSize*winSize/2+(i+1)*winSizeComb1;
		} else {
			for(int i=0; i<traceMaxSize; i++) {
				traceIndexContainer[i]=(i+1)*i*winSize/2+(i+1)*winSizeComb1;	
				

			}
		}

		// verified: a[] is set correctly.

		mat = initSumOfDistances(nse1, nse2, winSize, winSizeComb1, ca1, ca2);


		
//		try {
//			Matrix m2 = new Matrix(mat).copy();
//			JPanel panel = GuiWrapper.getScaleableMatrixPanel(m2);
//			JFrame frame = new JFrame();
//			frame.addWindowListener(new WindowAdapter(){
//				public void windowClosing(WindowEvent e){
//					JFrame f = (JFrame) e.getSource();
//					f.setVisible(false);
//					f.dispose();
//				}				
//			});
//						
//			
//			frame.getContentPane().add(panel);
//
//			frame.pack();
//			frame.setVisible(true);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

		
		// Set the distance matrix
		//afpChain.setDistanceMatrix(new Matrix(mat.clone()));


		//		
		//			   double rmsdThr = params.getRmsdThr();
		//			   StringBuffer buf = new StringBuffer("  ");  
		//			   for(int i=0; i<nse2; i++) 
		//			      buf.append(String.format("%c", i%10==0?(i%100)/10+48:32));
		//			   buf.append("\n");  
		//			   for(int i=0; i<nse1; i++) {
		//			      buf.append(String.format("%c ", i%10==0?(i%100)/10+48:32));
		//			      for(int j=0; j<nse2; j++) 
		//			         buf.append(String.format("%c", (mat[i][j])<rmsdThr?'+':'X'));
		//			      //printf("%c", ((int)*(mat[i]+j)/40)>9?'*':((int)*(mat[i]+j)/40)+48);
		//			      buf.append("\n");
		//			   }
		//			   buf.append("\n");
		//
		//			   System.out.println(buf.toString());
		//			


		return afpChain;
	}

	/**
	 * Evaluates the distance between two atoms
	 * Several scoring functions are implemented and can be changed by calling 
	 * {@link CeParameters#setScoringStrategy(Integer) setScoringStrategy()}
	 * on {@link CeParameters parameter} object this CECalculator was created with.
	 * <p>
	 * Scoring Strategies:<dl>
	 * <dt>DEFAULT_SCORING_STRATEGY</dt>
	 * <dd>Strategy of the original CE publication; CA-CA distance</dd>
	 * 
	 * <dt>SIDE_CHAIN_SCORING</dt>
	 * <dd>CB-CB distance. This performs better for sheets and helices than CA.</dd>
	 * 
	 * <dt>SIDE_CHAIN_ANGLE_SCORING</dt>
	 * <dd>Use the dot product (eg the cosine) of the two CA-CB vectors.</dd>
	 * 
	 * <dt>CA_AND_SIDE_CHAIN_ANGLE_SCORING</dt>
	 * <dd>Equivalent to DEFAULT_SCORING_STRATEGY + SIDE_CHAIN_ANGLE_SCORING</dd>
	 * </dl>
	 * 
	 *  <dt>SEQUENCE_CONSERVATION</dt>
	 * <dd>A mix between the DEFAULT_SCORING_STRATEGY and a scoring function that favors the alignment of sequence conserved positions in the alignment</dd>
	 * </dl>
	 * 
	 *   
	 * 
	 * @param ca1 The CA of the first residue
	 * @param ca2 The CA of the second residue
	 * @return The distance between the two fragments, according to the selected
	 * scoring strategy. Lower distances are better alignments.
	 * @throws StructureException
	 */
	private double getDistanceWithSidechain(Atom ca1, Atom ca2) throws StructureException {
		
		if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.CA_SCORING) {
			
			return Calc.getDistance(ca1,ca2);

		}

		double dist;
		Group g1 = ca1.getGroup();
		Atom cb1 = null;
		if ( g1.hasAtom(StructureTools.CB_ATOM_NAME)) {
			cb1 = g1.getAtom(StructureTools.CB_ATOM_NAME);
		}
		//
		Group g2 = ca2.getGroup();
		Atom cb2 = null;
		if ( g2.hasAtom(StructureTools.CB_ATOM_NAME)) {
			cb2 = g2.getAtom(StructureTools.CB_ATOM_NAME);
		}


		if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.SIDE_CHAIN_SCORING) {


			// here we are using side chain orientation for scoring...


			// score type 1    consider side chain distances   
			if ( cb1 != null && cb2 != null) {
				// CB distance
				dist = Calc.getDistance(cb1,cb2);
				//dist = dist / 2.;
			} else {
				dist = Calc.getDistance(ca1,ca2);
			}

			return dist;
		}

		else if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.SIDE_CHAIN_ANGLE_SCORING){

			// score type 2 add angle info


			if ( cb1 != null && cb2 != null) {
				// If the CA were overlaid, what is the distance between the CB?
				// Recall c^2 = a^2 + b^2 -2ab*cos(theta), so this is a function of angle
				Atom c1 = Calc.subtract(cb1, ca1);
				Atom c2 = Calc.subtract(cb2, ca2);
				Atom newA = Calc.subtract(c2, c1);
				dist = Calc.amount(newA); 
			}  else {
				//dist += Calc.getDistance(ca1,ca2);
				dist = 0;
			}

			return dist;

		}
		else if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.CA_AND_SIDE_CHAIN_ANGLE_SCORING){

			// score type 3
			// CA distance + cos(angle)
			dist = 0;
			if ( cb1 != null && cb2 != null) {
				Atom cacb1 = Calc.subtract(cb1, ca1);
				Atom cacb2 = Calc.subtract(cb2, ca2);
				Atom newA = Calc.subtract(cacb2, cacb1);
				//System.out.format("CACB 1: %s\nCACB 2: %s\ndiff: %s\nd: %f\n",cacb1.toString(),cacb2.toString(),newA.toString(),Calc.amount(newA));
				dist += Calc.amount(newA);
			}
			dist += Calc.getDistance(ca1,ca2);

			return dist;

		} else if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.SEQUENCE_CONSERVATION){
			if ( cb1 != null && cb2 != null) {
				// CB distance
				dist = Calc.getDistance(cb1,cb2);
				//dist = dist / 2.;
			} else {
				dist = Calc.getDistance(ca1,ca2);
			}
			return dist;
			
			
		}
		else {
			// unsupported scoring scheme
			return Calc.getDistance(ca1,ca2);
		}
	}

	/** build up intramolecular distance matrix dist1 & dist2
	 * 
	 * @param ca
	 * @param nse
	 * @return
	 * @throws StructureException
	 */
	private double[][] initIntraDistmatrix(Atom[] ca, int nse) throws StructureException
	{


		double[][] intraDist = new double[nse][nse];

		// 
		for(int ise1=0; ise1<nse; ise1++)  {

			for(int ise2=0; ise2<nse; ise2++)  {				
				intraDist[ise1][ise2] = getDistanceWithSidechain(ca[ise1], ca[ise2]);

			}
		}
		return intraDist;
	}


	public double[][] initSumOfDistances(int nse1, int nse2, int winSize, int  winSizeComb1, Atom[] ca1, Atom[] ca2) {

		double d;

		double[][] mat   = new double[nse1][nse2];

		// init the initial mat[] array.
		// at this stage mat contains the sum of the distances of fragments of the matrices dist1, dist
		for(int ise1=0; ise1<nse1; ise1++) {
			for(int ise2=0; ise2<nse2; ise2++) {

				mat[ise1][ise2]=-1.0;

				if(ise1>nse1-winSize || ise2>nse2-winSize) continue;

				d=0.0; 
				// this sums up over the distances of the fragments
				for(int is1=0; is1<winSize-2; is1++)
					for(int is2=is1+2; is2<winSize; is2++) {
						//System.out.println("pos1 :" +  (ise1+is1) + " " + (ise1+is2) +  " " + (ise2+is1) + " " + (ise2+is2));
						// is this abs or floor? check!
						d+=Math.abs(dist1[ise1+is1][ise1+is2]-dist2[ise2+is1][ise2+is2]);							
					}					
				mat[ise1][ise2]=d/winSizeComb1;					

				//System.out.println("mat ["+ise1+"]["+ise2+"]="+mat[ise1][ise2]);
			}

		}

		// verified: mat[][] probably ok.

		return mat;
	}






	@SuppressWarnings("unused")
	public void traceFragmentMatrix( AFPChain afpChain,
			Atom[] ca1, Atom[] ca2) {

		double rmsdThr = params.getRmsdThr();


		double oldBestTraceScore=10000.0;
		bestTraceScore = 100.0;
		nBestTrace=0;
		int nBestTrace0 = 0;
		int winSize = params.getWinSize();
		int winSizeComb1=(winSize-1)*(winSize-2)/2;
		boolean distAll = false;

		int winSizeComb2=distAll?winSize*winSize:winSize;
		double rmsdThrJoin = params.getRmsdThrJoin();

		double z0;		

		//double bestTraceZScore=-1.0;

		int nse1 = ca1.length;
		int nse2 = ca2.length;

		//System.out.println("nse1 :" +nse1 + " nse2: " + nse2);

		int traceMaxSize=nse1<nse2?nse1:nse2;

		bestTrace1 = new int [traceMaxSize];
		bestTrace2 = new int [traceMaxSize];
		trace1     = new int [traceMaxSize]; 
		trace2     = new int [traceMaxSize];

		int[] traceIndex     = new int [traceMaxSize];
		int[] traceIterLevel = new int [traceMaxSize];

		int ise11;
		int ise12;
		int ise21;
		int ise22;

		int ise1;
		int ise2;

		int gapMax=params.getMaxGapSize();

		int iterDepth;
		if ( gapMax > 0){
			iterDepth 	=gapMax*2+1;
		} else {
			iterDepth = traceMaxSize;
		}
		double[][] traceScore = new double[traceMaxSize][iterDepth];

		nTraces =0;
		long tracesLimit=(long)5e7;
		double score =-1;
		double score0 = -1;
		double score1 = -1 ;
		double score2 = -1;

		int mse1;
		int mse2;
		int jgap;
		int jdir;
		int jse1=0;
		int jse2=0;

		int bestTracesMax=30;
		bestTraces1 = new int[bestTracesMax][traceMaxSize];
		bestTraces2 = new int[bestTracesMax][ traceMaxSize];
		bestTracesN=new int [bestTracesMax];
		bestTracesScores = new double [bestTracesMax];
		for(int it=0; it<bestTracesMax; it++) {
			bestTracesN[it]=0;
			bestTracesScores[it]=100;
		}

		nBestTraces=0;
		int newBestTrace=0;

		double traceTotalScore=0;
		double traceScoreMax =0;
		double userRMSDMax = params.getMaxOptRMSD();
		int kse1;
		int kse2;

		iterLoop:
			for(int iter=0; iter<nIter; iter++) {

				if(iter>2) {
					if(oldBestTraceScore<=bestTraceScore) break;
				}
				oldBestTraceScore=bestTraceScore;

				if(iter==1) {
					z0=zStrAlign(winSize, nBestTrace, bestTraceScore, 
							bestTrace1[nBestTrace]+winSize-bestTrace1[0]+
							bestTrace2[nBestTrace]+winSize-bestTrace2[0]-
							nBestTrace*2*winSize);
					if(z0<zThr) break;
					nBestTrace0=nBestTrace;
					nBestTrace=0; 
					bestTraceScore=100.0;

					nTraces=0;
				}


				if(iter==0) {
					ise11=0; ise12=nse1;
					ise21=0; ise22=nse2;

				}
				else {
					if(iter==1) {
						ise11=bestTrace1[0]; ise12=bestTrace1[0]+1;
						ise21=bestTrace2[0]; ise22=bestTrace2[0]+1;
					}
					else {
						ise11=bestTrace1[0]-1; ise12=bestTrace1[0]+2;
						ise21=bestTrace2[0]-1; ise22=bestTrace2[0]+2;
					}
					if(ise11<0) ise11=0;
					if(ise12>nse1) ise12=nse1;
					if(ise21<0) ise21=0;
					if(ise22>nse2) ise22=nse2;				
				}

				//System.out.println("ise1Loop: " + ise11 + " " + ise12 + " " + ise21 + " " + ise22);
				ise1Loop:
					for(int ise1_=ise11; ise1_<ise12; ise1_++) {
						ise2Loop:
							for(int ise2_=ise21; ise2_<ise22; ise2_++) {

								ise1=ise1_;
								ise2=ise2_;
								if(iter>1 && ise1==ise11+1 && ise2==ise21+1) continue ise2Loop;

								//if(ise2==ise21) 	System.out.println(String.format("(%d, %d)",ise1, nTraces));


								if(iter==0 && (ise1>nse1-winSize*(nBestTrace-1) || 
										ise2>nse2-winSize*(nBestTrace-1))) continue ise2Loop;

								if(mat[ise1][ise2]<0.0) continue ise2Loop;
								if(mat[ise1][ise2]>rmsdThr) continue ise2Loop;
								if (mat[ise1][ise2]>userRMSDMax) continue ise2Loop;
								nTrace=0;
								trace1[nTrace]=ise1; 
								trace2[nTrace]=ise2;
								traceIndex[nTrace]=0;
								traceIterLevel[nTrace]=0;

								score0=mat[ise1][ise2];


								nTrace++;
								boolean isTraceUp=true;
								int traceIndex_=0;

								traceLoop:
									while(nTrace>0) {

										kse1=trace1[nTrace-1]+winSize;
										kse2=trace2[nTrace-1]+winSize;

										//System.out.println("isTraceUp " + isTraceUp + " " + nTrace + " " + kse1 + " " + kse2);

										while(true) {
											if(kse1>nse1-winSize-1) break;
											if(kse2>nse2-winSize-1) break;
											if(mat[kse1][kse2]>=0.0) break;
											kse1++;
											kse2++;
										}


										traceIndex_=-1; 

										if(isTraceUp) {

											int nBestExtTrace=nTrace;
											double bestExtScore=100.0;


											// extension of the alignment path
											// condition 4, 5
											itLoop:
												for(int it=0; it<iterDepth; it++) {

													jgap=(it+1)/2;
													jdir=(it+1)%2;

													if(jdir==0) {
														mse1=kse1+jgap;
														mse2=kse2;
													}
													else {
														mse1=kse1;
														mse2=kse2+jgap;
													}

													if(mse1>nse1-winSize-1) continue itLoop;
													if(mse2>nse2-winSize-1) continue itLoop;

													if(mat[mse1][mse2]<0.0)     continue itLoop;
													if(mat[mse1][mse2]>rmsdThr) continue itLoop;
													if(mat[mse1][mse2]>userRMSDMax) continue itLoop;
													
													nTraces++;
													if(nTraces>tracesLimit) {

														return;
													}

													score=0.0;

													//													if(!distAll) {
														//System.out.println("getting score " + mse1 + " " + mse2 + " " + winSize + " " + jgap + " " + jdir + " " + it + " " + kse1 + " " + kse2);
														score = getScoreFromDistanceMatrices(mse1,mse2,winSize);
														//System.out.println("got score: " + score);
														score1=score/(nTrace*winSize);

														//													} else {
															//														// all dist
															//														for(int itrace=0; itrace<nTrace; itrace++) {
														//															for(int is1=0; is1<winSize; is1++)
														//																for(int is2=0; is2<winSize; is2++)
														//																	score+=Math.abs(dist1[trace1[itrace]+is1][mse1+is2]-
														//																			dist2[trace2[itrace]+is1][mse2+is2]);
														//														}
														//														score1=score/(nTrace*winSize*winSize);
														//													}


														//System.out.println("up: " + nTrace + " "  + score + " " + score0 + " " + score1 + " " + winSize + " " + traceIndex_ + " " + it + " ");
														if(score1>rmsdThrJoin) 
															continue itLoop;
														if(score1>userRMSDMax)
														   continue itLoop;
														score2=score1;

														// this just got checked, no need to check again..
														//if(score2>rmsdThrJoin) 
														//	continue itLoop;

														if(nTrace>nBestExtTrace || (nTrace==nBestExtTrace &&
																score2<bestExtScore)) {
															//System.out.println("setting traceindex to " + it + " " + score2);
															bestExtScore=score2;
															nBestExtTrace=nTrace;
															traceIndex_=it;
															traceScore[nTrace-1][traceIndex_]=score1;
														}

												}
										}

										if(traceIndex_!=-1) {
											jgap=(traceIndex_+1)/2;
											jdir=(traceIndex_+1)%2;
											if(jdir==0) {
												jse1=kse1+jgap;
												jse2=kse2;
											}
											else {
												jse1=kse1;
												jse2=kse2+jgap;
											}

											if(iter==0){

												score1=(traceScore[nTrace-1][traceIndex_]*winSizeComb2*nTrace+
														mat[jse1][jse2]*winSizeComb1)/(winSizeComb2*nTrace+
																winSizeComb1);

												score2 = getScore2(jse1, jse2, traceScore, traceIndex_, traceIndex, winSizeComb1, winSizeComb2, score0, score1);

												if(score2>rmsdThrJoin) 
													traceIndex_=-1;
												else if ( score2 > userRMSDMax) 
												   traceIndex_=-1;												
												else {
													traceScore[nTrace-1][traceIndex_]=score2;

													traceTotalScore=score2;
												}

											}
											else {
												if(traceScoreMax>rmsdThrJoin && nBestTrace>=nBestTrace0) 
													traceIndex_=-1;								
												traceTotalScore=traceScoreMax;
											}
										}

										//System.out.println("middle: " + nTrace + " " + score + " " + score0 + " " + score1 + "  " + score2  + " " + traceIndex_);

										if(traceIndex_==-1) {
											//System.out.println("continue traceLoop " + nTrace);
											//if(iterLevel==1) break;
											nTrace--;
											isTraceUp=false;
											continue traceLoop;
										}
										else {
											traceIterLevel[nTrace-1]++;
											trace1[nTrace]=jse1;
											trace2[nTrace]=jse2;
											traceIndex[nTrace]=traceIndex_;
											traceIterLevel[nTrace]=0;
											nTrace++;
											isTraceUp=true;

											if(nTrace>nBestTrace || 
													(nTrace==nBestTrace  && 
															bestTraceScore>traceTotalScore)) {

												for(int itrace=0; itrace<nTrace; itrace++) {
													bestTrace1[itrace]=trace1[itrace];
													bestTrace2[itrace]=trace2[itrace];
												}
												bestTraceScore=traceTotalScore;
												nBestTrace=nTrace;
											}

											if(iter==0) {
												//System.out.println("doing iter0 " + newBestTrace + " " + traceTotalScore + " " + bestTracesMax);
												newBestTrace = doIter0(newBestTrace,traceTotalScore, bestTracesMax);


											}
										}
									}
							}
					}

				if ( isPrint) {
					System.out.println("fragment length: " + params.getWinSize());
					System.out.println("ntraces : " + nTraces );
				}



			}

//		try {
//			Matrix m2 = new Matrix(traceScore).copy();
//			JPanel panel = GuiWrapper.getScaleableMatrixPanel(m2);
//			JFrame frame = new JFrame();
//			frame.addWindowListener(new WindowAdapter(){
//				public void windowClosing(WindowEvent e){
//					JFrame f = (JFrame) e.getSource();
//					f.setVisible(false);
//					f.dispose();
//				}				
//			});
//						
//			
//			frame.getContentPane().add(panel);
//
//			frame.pack();
//			frame.setVisible(true);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}


		if ( params.isShowAFPRanges()){
			System.out.println("fragment length: " + params.getWinSize());
			System.out.println("ntraces : " + nTraces );         

		}

	}

	protected double getScore2(int jse1, int jse2, double[][] traceScore, int traceIndex_,int[] traceIndex,int winSizeComb1, int winSizeComb2, double score0, double score1 ) {



		/*double score2=
			((nTrace>1?traceScore[nTrace-2][traceIndex[nTrace-1]]:score0)
		 *a[nTrace-1]+score1*(a[nTrace]-a[nTrace-1]))/a[nTrace];
		 */
		double val = 0;
		if ( nTrace>1)
			val =traceScore[nTrace-2][traceIndex[nTrace-1]];
		else
			val = score0;

		double score2 =  (val * traceIndexContainer[nTrace-1]+score1*(traceIndexContainer[nTrace]-traceIndexContainer[nTrace-1]))/traceIndexContainer[nTrace];

		//System.out.println("check: score0 " + score0 + " score 1 " + score1 + " sc2: " + score2 + " val: " + val + " nTrace:" + nTrace+ " " +  traceIndexContainer[nTrace-1] + " " +  traceIndexContainer[nTrace-1] + " " + traceIndexContainer[nTrace] );

		return score2;


	}

	protected int doIter0(int newBestTrace, double traceTotalScore, double bestTracesMax) {


		// only do the whole method if this criteria is fulfilled...
		if(nTrace>bestTracesN[newBestTrace] ||
				(nTrace==bestTracesN[newBestTrace] && 
						bestTracesScores[newBestTrace]>traceTotalScore)) {


			for(int itrace=0; itrace<nTrace; itrace++) {
				bestTraces1[newBestTrace][itrace]=trace1[itrace];
				bestTraces2[newBestTrace][itrace]=trace2[itrace];
				bestTracesN[newBestTrace]=nTrace;
				bestTracesScores[newBestTrace]=traceTotalScore;
				//System.out.println("bestTracesScrores ["+newBestTrace+"]=" +traceTotalScore);
			}

			if(nTrace>nBestTrace) nBestTrace=nTrace;

			if(nBestTraces<bestTracesMax) {
				nBestTraces++;
				newBestTrace++;
			}

			if(nBestTraces==bestTracesMax) {
				//System.out.println("nBestTraces == bestTracesmax " + nBestTraces);
				newBestTrace=0;
				double scoreTmp=bestTracesScores[0];
				int nTraceTmp=bestTracesN[0];
				for(int ir=1; ir<nBestTraces; ir++) {
					if(bestTracesN[ir]<nTraceTmp || 
							(bestTracesN[ir]==nTraceTmp && 
									scoreTmp<bestTracesScores[ir])) {
						nTraceTmp=bestTracesN[ir];
						scoreTmp=bestTracesScores[ir];
						newBestTrace=ir;
						//System.out.println("setting new bestTracesScore to " + ir + " " + scoreTmp);
					}
				}
			}
		}

		//System.out.println("iter0 : " + newBestTrace + " " + bestTracesN[newBestTrace] + " " + traceTotalScore + " " + nTrace);



		/*
z=zStrAlign(winSize, nTrace, traceTotalScore, 
trace1[nTrace-1]-trace1[0]+trace2[nTrace-1]-trace2[0]-
2*(nTrace-1)*winSize);
if(z>bestTraceZScore) {
for(int itrace=0; itrace<nTrace; itrace++) {
bestTrace1[itrace]=trace1[itrace];
bestTrace2[itrace]=trace2[itrace];
}
bestTraceZScore=z;
bestTraceScore=*(traceScore[nTrace-2]+traceIndex_);
nBestTrace=nTrace;
}
		 */
		return newBestTrace;

	}


	protected double getScoreFromDistanceMatrices(int mse1, int mse2,int winSize) {

		double score = 0;
		// (winSize) "best" dist

		// reduce sign. values to C code.. 6 digits..

		for(int itrace=0; itrace<nTrace; itrace++) {
			score+=  Math.abs(dist1[trace1[itrace]][mse1]-
					dist2[trace2[itrace]][mse2]);

			score+=  Math.abs(dist1[trace1[itrace]+winSize-1]
			                        [mse1+winSize-1]-
			                        dist2[trace2[itrace]+winSize-1][mse2+winSize-1]);

			for(int id=1; id<winSize-1; id++) 
				score+=  Math.abs(dist1[trace1[itrace]+id][mse1+winSize-1-id]-
						dist2[trace2[itrace]+id][mse2+winSize-1-id]);

		}
		
		return score;
	}

	public void nextStep( AFPChain afpChain,
			Atom[] ca1, Atom[] ca2) throws StructureException{


		if(nBestTrace>0) {
			checkBestTraces(afpChain,ca1,ca2);
		} else {
			noBestTrace();
		}

		convertAfpChain(afpChain, ca1, ca2);
		AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2);
	}



	// this part is modified from the original CeCalculator
	@SuppressWarnings("unused")
	private void checkBestTraces( AFPChain afpChain,
			Atom[] ca1, Atom[] ca2) throws StructureException{

		z=0.0;


		int nGaps;
		int winSize = params.getWinSize();
		int nse1 = ca1.length;
		int nse2 = ca2.length;
		int traceMaxSize=nse1<nse2?nse1:nse2;
		int idir;


		align_se1=new int [nse1+nse2];
		align_se2=new int [nse1+nse2];
		lcmp = 0;

		// we now support alignment using any particular atoms..

		Atom[] strBuf1 = new Atom[traceMaxSize];
		Atom[] strBuf2 = new Atom[traceMaxSize];

		double rmsdNew;



		// removing some loops that are run in orig CE
		// and which did not do anything
		if ( debug ){			
			checkPrintRmsdNew(traceMaxSize, winSize, ca1, ca2);						
		}

		double rmsd=100.0;

		int iBestTrace=0;

		for(int ir=0; ir<nBestTraces; ir++) {
			if(bestTracesN[ir]!=nBestTrace) continue;

			rmsdNew = getRMSDForBestTrace(ir, strBuf1, strBuf2, bestTracesN,bestTraces1, bestTrace2,winSize,ca1,ca2);
			if ( isPrint)
				System.out.println(String.format("%d %d %d %.2f", ir, bestTracesN[ir], nBestTrace, rmsdNew));

			if(rmsd>rmsdNew) {
				iBestTrace=ir;
				rmsd=rmsdNew;
				//System.out.println(" iBestTrace:" + iBestTrace + " new rmsd = " + rmsd);
			}
		}
		for(int it=0; it<bestTracesN[iBestTrace]; it++) {
			bestTrace1[it]=bestTraces1[iBestTrace][it];
			bestTrace2[it]=bestTraces2[iBestTrace][it];
		}

		//System.out.println("iBestTrace: "+iBestTrace+" = bestTracesScores " + bestTracesScores[iBestTrace]);

		nBestTrace=bestTracesN[iBestTrace];

		bestTraceScore=bestTracesScores[iBestTrace];


		//printf("\nOptimizing gaps...\n");

		int[] traceLen=new int [traceMaxSize];
		bestTraceLen=new int [traceMaxSize];


		int strLen=0;

		int jt;
		strLen=0;
		nGaps=0;    
		nTrace=nBestTrace;  

		for(jt=0; jt<nBestTrace; jt++) {
			trace1[jt]=bestTrace1[jt];
			trace2[jt]=bestTrace2[jt];
			traceLen[jt]=winSize;

			if(jt<nBestTrace-1) {
				nGaps+=bestTrace1[jt+1]-bestTrace1[jt]-winSize+
				bestTrace2[jt+1]-bestTrace2[jt]-winSize;
			}
		}    
		nBestTrace=0;
		for(int it=0; it<nTrace; ) {
			int cSize=traceLen[it];
			for(jt=it+1; jt<nTrace; jt++) {
				if(trace1[jt]-trace1[jt-1]-traceLen[jt-1]!=0 ||
						trace2[jt]-trace2[jt-1]-traceLen[jt-1]!=0) break;
				cSize+=traceLen[jt]; 
			}
			bestTrace1[nBestTrace]=trace1[it];
			bestTrace2[nBestTrace]=trace2[it];
			bestTraceLen[nBestTrace]=cSize;
			nBestTrace++;
			strLen+=cSize;
			it=jt;
		}


		int is=0;
		for(jt=0; jt<nBestTrace; jt++) {
			for(int i=0; i<bestTraceLen[jt]; i++) {
				setStrBuf(strBuf1,is+i,ca1,bestTrace1[jt]+i );
				setStrBuf(strBuf2,is+i,ca2,bestTrace2[jt]+i);
			}
			is+=bestTraceLen[jt];
		}
		//sup_str(strBuf1, strBuf2, strLen, d_);

		rmsd=calc_rmsd(strBuf1, strBuf2, strLen,true);

		if ( isPrint)
			System.out.println("got first rmsd: " + rmsd);
		boolean isCopied=false;

		outer_loop:
			for(int it=1; it<nBestTrace; it++) {			

				/* not needed...
			int igap;
			if(bestTrace1[it]-bestTrace1[it-1]-bestTraceLen[it-1]>0) igap=0;
			if(bestTrace2[it]-bestTrace2[it-1]-bestTraceLen[it-1]>0) igap=1;
				 */


				boolean wasBest=false;
				main_loop:
					for(idir=-1; idir<=1; idir+=2) {
						if(wasBest) break;

						inner_loop:
							for(int idep=1; idep<=winSize/2; idep++) {

								// isCopied indicates that bestTrace has changed and needs to be re-copied
								if(!isCopied)
									for(jt=0; jt<nBestTrace; jt++) {
										trace1[jt]=bestTrace1[jt];
										trace2[jt]=bestTrace2[jt];
										traceLen[jt]=bestTraceLen[jt];
									}
								isCopied=false;

								// Move an atom from the previous trace to the current on, or vice versa
								traceLen[it-1]+=idir;
								traceLen[it]-=idir;
								trace1[it]+=idir;
								trace2[it]+=idir;

								// Copy atoms from the current trace into strBuf
								is=0;
								for(jt=0; jt<nBestTrace; jt++) {
									for(int i=0; i<traceLen[jt]; i++) {
										if(ca1[trace1[jt]+i].getX()>1e10 || ca2[trace2[jt]+i].getX()>1e10) 
											continue main_loop;
										strBuf1[is+i]=ca1[trace1[jt]+i];
										strBuf2[is+i]=ca2[trace2[jt]+i];
									}
									is+=traceLen[jt];
								}
								// Check new RMSD
								//sup_str(strBuf1, strBuf2, strLen, d_);
								rmsdNew=calc_rmsd(strBuf1, strBuf2, strLen, true);
								//System.out.println(String.format("step %d %d %d %.2f old: %.2f", it, idir, idep, rmsdNew, rmsd));
								
								// Update best trace if RMSD improved
								if(rmsdNew<rmsd) {

									for(jt=0; jt<nBestTrace; jt++) {
										bestTrace1[jt]  = trace1[jt];
										bestTrace2[jt]  = trace2[jt];
										bestTraceLen[jt]= traceLen[jt];
									}
									isCopied=true;
									wasBest=true;
									rmsd=rmsdNew;
									continue inner_loop;
								}
								// AP
								//bad_ca: break;
								continue main_loop;
							}
					}
			}
		rmsdNew=calc_rmsd(strBuf1, strBuf2, strLen,true);
		if ( isPrint)
			System.out.println("rmsdNew: " + rmsdNew + " rmsd " + rmsd);
		afpChain.setTotalRmsdIni(rmsdNew);
		afpChain.setTotalLenIni(strBuf1.length);


		nAtom = strLen;

		//System.out.println("zStrAlign: " + winSize + " strLen " + strLen  + " s/w " + (strLen/winSize) + " " + bestTraceScore + " " + nGaps);
		z=zStrAlign(winSize, strLen/winSize, bestTraceScore, nGaps);

		if(params.isShowAFPRanges()) {
			System.out.println("win size: " + winSize + " strLen/winSize: " + strLen/winSize + " best trace score: " + String.format("%.2f",bestTraceScore) + " nr gaps: " + nGaps + " nr residues: " + nAtom);

			System.out.println(String.format("size=%d rmsd=%.2f z=%.1f gaps=%d(%.1f%%) comb=%d", 
					nAtom, rmsd, z, nGaps, nGaps*100.0/nAtom,
					nTraces)); 

			System.out.println("Best Trace, before optimization");
			for(int k=0; k<nBestTrace; k++)
				System.out.println(String.format("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,
						bestTraceLen[k]));

		}

		// start to convert CE internal datastructure to generic AFPChain one...
		List<AFP> afpSet = new ArrayList<AFP>();
		for (int afp=0;afp<nBestTrace;afp++){
			// fill in data from nBestTrace into AFP

			AFP afpI = new AFP();

			afpI.setFragLen(bestTraceLen[afp]);
			afpI.setP1(bestTrace1[afp]+1);
			afpI.setP2(bestTrace2[afp]+1);

			afpSet.add( afpI);
		}

		afpChain.setAfpSet(afpSet);




		//System.out.println("z:"+z + " zThr" + zThr+ " bestTraceScore " + bestTraceScore + " " + nGaps );
		if(z>=zThr) {
			nGaps = optimizeSuperposition(afpChain,nse1, nse2, strLen, rmsd, ca1, ca2,nGaps,strBuf1,strBuf2);
			//	      if(isPrint) {
			//		/*
			//		FILE *f=fopen("homologies", "a");
			//		fprintf(f, "%s(%d) %s(%d) %3d %4.1f %4.1f %d(%d) ", 
			//			name1, nse1, name2, nse2, nAtom, rmsd, z, 
			//			nGaps, nGaps*100/nAtom);
			//		for(int k=0; k<nBestTrace; k++)
			//		  fprintf(f, "(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1, 
			//			  bestTraceLen[k]);
			//		fprintf(f, "\n");
			//		fclose(f);
			//		*/
			//	      }
		}
		else {
			int lali_x_ = 0;
			for(int k=0; k<nBestTrace; k++) {
				for(int l=0; l<bestTraceLen[k]; l++) {
					align_se1[lcmp+l]=bestTrace1[k]+l;
					align_se2[lcmp+l]=bestTrace2[k]+l;
				}
				lali_x_+=bestTraceLen[k];
				if(k<nBestTrace-1) {
					if(bestTrace1[k]+bestTraceLen[k]!=bestTrace1[k+1])
						for(int l=bestTrace1[k]+bestTraceLen[k]; l<bestTrace1[k+1]; l++) {
							align_se1[lcmp]=l;
							align_se2[lcmp]=-1;
							lcmp++;
						}
					if(bestTrace2[k]+bestTraceLen[k]!=bestTrace2[k+1])
						for(int l=bestTrace2[k]+bestTraceLen[k]; l<bestTrace2[k+1]; l++) {
							align_se1[lcmp]=-1;
							align_se2[lcmp]=l;
							lcmp++;
						}
				}
			}
			nAtom=lali_x_;
		}

		timeEnd = System.currentTimeMillis();
		long time_q=(timeEnd-timeStart);

		double gapsP = ( nGaps*100.0/nAtom) ;
		if(isPrint) {
			String msg = String.format("Alignment length = %d Rmsd = %.2fA Z-Score = %.1f Gaps = %d(%.1f%%)",nAtom,rmsd,z,nGaps, gapsP);
			System.out.println(msg + " CPU = " + time_q);
		}

		//      if ( params.isShowAFPRanges()){

		// this is actually the final alignment...
		//         System.out.println("Best Trace: (index1,index2,len)");
		//         for(int k=0; k<nBestTrace; k++)
		//            System.out.println(
		//                  String.format("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1, bestTraceLen[k]));
		//
		//
		//
		//      }

		afpChain.setCalculationTime(time_q);
		afpChain.setGapLen(nGaps);

		int[] optLen = new int[]{nAtom};
		afpChain.setOptLen(optLen);
		afpChain.setOptLength(nAtom);		
		afpChain.setAlnLength(lcmp);

		afpChain.setProbability(z);


	}

	/** set the Atoms for a particular residue position.
	 * Requires that atom.getParent returns the correct group!
	 * take care during cloning of atoms. Best to use StructureTools.cloneCaAtoms();
	 * 
	 * @param strBuf
	 * @param i
	 * @param ca
	 * @param j
	 */
	private void setStrBuf(Atom[] strBuf, int i, Atom[] ca, int j) {
		// TODO Auto-generated method stub
		//TODO
		Group parent = ca[j].getGroup();
		int pos = 0;
		String atomName = ca[j].getName();

		Atom a = null;
		
			a= parent.getAtom(atomName);
			if ( a != null){
				strBuf[i]=a;
			} else {
				// 	probably a GLY and no CB was found...
				//e.printStackTrace();
			}
		strBuf[i+pos] = a;
		pos++;



	}

	// TODO:  consider all requested Atoms?
	private double getRMSDForBestTrace(int ir, Atom[] strBuf1, Atom[] strBuf2,
			int[] bestTracesN2, int[][] bestTraces12, int[] bestTrace22, 
			int winSize,Atom[] ca1, Atom[] ca2 ) throws StructureException {
		int is=0;
		for(int jt=0; jt<bestTracesN[ir]; jt++) {
			for(int i=0; i<winSize; i++) {

				setStrBuf(strBuf1, is+i, ca1, bestTraces1[ir][jt]+i);
				setStrBuf(strBuf2, is+i, ca2, bestTraces2[ir][jt]+i);
			}
			is+=winSize;
		}
		//sup_str(strBuf1, strBuf2, bestTracesN[ir]*winSize, d_);
		double rmsdNew=calc_rmsd(strBuf1, strBuf2, bestTracesN[ir]*winSize, true);
		return rmsdNew;

	}






	/** calc initial RMSD for bestTrace1 in debug only
	 * 
	 */
	private void checkPrintRmsdNew(int traceMaxSize, int winSize, Atom[] ca1, Atom[] ca2) throws StructureException{

		int is = 0;
		// calc initial RMSD for bestTrace1
		Atom[] strBuf1 = new Atom[traceMaxSize];
		Atom[] strBuf2 = new Atom[traceMaxSize];
		for(int jt=0; jt<nBestTrace; jt++) {
			for(int i=0; i<winSize; i++) {
				setStrBuf(strBuf1, is+i, ca1, bestTrace1[jt]+i );
				setStrBuf(strBuf2, is+i, ca2, bestTrace2[jt]+i );
			}
			is+=winSize;
		}

		//sup_str(strBuf1, strBuf2, nBestTrace*winSize, d_);
		double rmsdNew=calc_rmsd(strBuf1, strBuf2, nBestTrace*winSize, true);
		//afpChain.setTotalRmsdIni(rmsdNew);


		if ( isPrint){
			System.out.println("rmsdNew after trace: " +rmsdNew);

			for(int k=0; k<nBestTrace; k++)
				System.out.println(String.format("(%d,%d,%d) ", bestTrace1[k]+1, bestTrace2[k]+1,8));
		}

		if ( isPrint){
			System.out.println("best traces: " + nBestTraces);
		}


	}






	private static char getOneLetter(Group g){

		if (g==null) return StructureTools.UNKNOWN_GROUP_LABEL;
		
		return StructureTools.get1LetterCode(g.getPDBName());
	}



	private int optimizeSuperposition(AFPChain afpChain, int nse1, int nse2, int strLen, double rmsd, Atom[] ca1, Atom[] ca2,int nGaps, 
			Atom[] strBuf1, Atom[] strBuf2 ) throws StructureException {

		//System.out.println("optimizing Superimposition...");

		//nAtom=strLen;
		// optimization on superposition
		Atom[] ca3=new Atom[nse2];


		double rmsdLen  = 0.0;

		// this flag tests if the RMSDLen has been assigned.
		// this is to enforce that the alignment ends up 
		// smaller than 95% of the original alignment.
		// +/- 
		boolean isRmsdLenAssigned=false;
		int nAtomPrev=-1;

		double oRmsdThr = params.getORmsdThr();
		
		double distanceIncrement = params.getDistanceIncrement();
		double maxUserRMSD = params.getMaxOptRMSD();
		nAtom=0;
		int counter = -1;
		
		int maxNrIterations = params.getMaxNrIterationsForOptimization();
		while((nAtom<strLen*0.95 || 
				(isRmsdLenAssigned && rmsd<rmsdLen*1.1 && nAtomPrev!=nAtom)) && ( counter< maxNrIterations)) {
		  
			counter++;
			if ( debug)
			   System.out.println("nAtom: " + nAtom + " " + nAtomPrev + " " + rmsdLen + " " + isRmsdLenAssigned + " strLen:" + strLen);
			nAtomPrev=nAtom;
			oRmsdThr += distanceIncrement;
			
			rot_mol(ca2, ca3, nse2, r,t);

			for(int ise1=0; ise1<nse1; ise1++) {
				for(int ise2=0; ise2<nse2; ise2++) {

					//mat[ise1][ise2]=-0.001;

					// this needs to be a parameter...


					double dist = getDistanceWithSidechain(ca1[ise1], ca3[ise2]);
					mat[ise1][ise2] = oRmsdThr - dist;
					
					//double distold = Calc.getDistance(ca1[ise1],ca3[ise2]);
					//double scoreOld  = oRmsdThr - distold ;
					//mat[ise1][ise2] = scoreOld;
					//mat[ise1][ise2] = oRmsdThr - Calc.getDistance(ca1[ise1],ca3[ise2]);

					//if ( counter == 0 &&  ise1 == ise2) {

					// System.out.println("mat[" + ise1 + "][" + ise2 + "] " + mat[ise1][ise2] + " scoreOld:" + scoreOld + " oRmsdThr: " + oRmsdThr +" dist: " + dist + " distold:" + distold );
					// }


				}
			}
			
			if ( params.getScoringStrategy() == CeParameters.ScoringStrategy.SEQUENCE_CONSERVATION){
				mat = updateMatrixWithSequenceConservation(mat,ca1,ca2, params);
			}
			
			mat = notifyMatrixListener(mat);
			
			double gapOpen = params.getGapOpen();
			double gapExtension = params.getGapExtension();
			
			double score = dpAlign( nse1, nse2, gapOpen , gapExtension , false, false);

			if (debug)
				System.out.println("iter: "+ counter + "  score:"  + score + " " + " nAtomPrev: " + nAtomPrev + " nAtom:" + nAtom + " oRmsdThr: " + oRmsdThr);

			afpChain.setAlignScore(score);


			nAtom=0; nGaps=0; 
			for(int ia=0; ia<lcmp; ia++) {
				if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {

					strBuf1[nAtom]=ca1[align_se1[ia]];
					strBuf2[nAtom]=ca2[align_se2[ia]];

					nAtom++;

				}
				else {
					nGaps++;
				}
			}



			if(nAtom<4) continue;

			//sup_str(strBuf1, strBuf2, nAtom, _d);
			// here we don't store the rotation matrix for the user!
			rmsd= calc_rmsd(strBuf1, strBuf2, nAtom,false);
			if ( isPrint )
				System.out.println("iter: " + counter + " nAtom " + nAtom + " rmsd: " + rmsd);
			//afpChain.setTotalRmsdOpt(rmsd);
			//System.out.println("rmsd: " + rmsd);
			
			if(!(nAtom<strLen*0.95) && (!isRmsdLenAssigned)) { 
				rmsdLen=rmsd;
				isRmsdLenAssigned=true;
			}	
			//System.out.println(String.format("nAtom %d %d rmsd %.1f", nAtom, nAtomPrev, rmsd));


			afpChain.setBlockRmsd(new double[]{rmsd});
			afpChain.setOptRmsd(new double[]{rmsd});
			afpChain.setTotalRmsdOpt(rmsd);
			afpChain.setChainRmsd(rmsd);
			
			if ( rmsd >= maxUserRMSD) {
	              break;
			}

		}



		//System.out.println("done optimizing");
		/*
		nAtom=0; nGaps=0; 
		for(int ia=0; ia<lcmp; ia++)
		if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
		if(ca1[align_se1[ia]].X<1e10 && ca2[align_se2[ia]].X<1e10) {
		strBuf1[nAtom]=ca1[align_se1[ia]];
		strBuf2[nAtom]=ca2[align_se2[ia]];
		nAtom++;
		}
		}
		else {
		nGaps++;
		}

		sup_str(strBuf1, strBuf2, nAtom, _d);
		rmsd=calc_rmsd(strBuf1, strBuf2, nAtom, _d);
		 */
		nBestTrace=0;
		boolean newBestTrace=true;
		for(int ia=0; ia<lcmp; ia++) {
			if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
				//System.out.println(" " +align_se1[ia] + " " + align_se2[ia]);

				if(newBestTrace) {						
					bestTrace1[nBestTrace]=align_se1[ia];
					bestTrace2[nBestTrace]=align_se2[ia];
					bestTraceLen[nBestTrace]=0;
					newBestTrace=false;
					nBestTrace++;
				}
				bestTraceLen[nBestTrace-1]++;

			}
			else {
				newBestTrace=true;
			}
		}

		return nGaps;

		// end of optimization on superposition

	}

	/** Modifies an alignment matrix by favoring the alignment of similar and identical amino acids and penalizing the alignment of unrelated ones.
	 * 
	 * @param max alignment matrix
	 * @param ca1 Atoms for protein 1
	 * @param ca2 Atoms for Protein 2
	 * @param params alignment parameters
	 * @return modified alignment matrix
	 */
	public static double[][] updateMatrixWithSequenceConservation(double[][] max, Atom[] ca1, Atom[] ca2, CeParameters params) {
		Matrix origM = new Matrix(max);
		
		SubstitutionMatrix<AminoAcidCompound> substMatrix =
			params.getSubstitutionMatrix();
		
		int internalScale = 1;
		if ( substMatrix instanceof ScaledSubstitutionMatrix) {
			ScaledSubstitutionMatrix scaledMatrix = (ScaledSubstitutionMatrix) substMatrix;
			internalScale = scaledMatrix.getScale();
		}

		
		AminoAcidCompoundSet set = AminoAcidCompoundSet.getAminoAcidCompoundSet();

		for (int i = 0 ; i < origM.getRowDimension() ; i++){
			for ( int j =0; j < origM.getColumnDimension() ; j ++ ) {
				double val = origM.get(i,j);
				Atom a1 = ca1[i];
				Atom a2 = ca2[j];

				AminoAcidCompound ac1 =
					set.getCompoundForString(a1.getGroup().getChemComp().getOne_letter_code());
				AminoAcidCompound ac2 =
					set.getCompoundForString(a2.getGroup().getChemComp().getOne_letter_code());
				
				
				if ( ac1 == null || ac2 == null)
					continue;
				
				short aaScore = substMatrix.getValue(ac1,ac2);
				
				double weightedScore = (aaScore / internalScale) * params.getSeqWeight();
				
				
				val += weightedScore;
				origM.set(i,j,val);

			}
		}
		max = origM.getArray();

		//SymmetryTools.showMatrix((Matrix)origM.clone(), "in optimizer "  + loopCount  );
		//SymmetryTools.showMatrix(origM, "iteration  matrix " + loopCount + " after");
		
		return max;
	}

	private double[][] notifyMatrixListener(double[][] mat2) {
		for (MatrixListener li : matrixListeners) {
			mat2 = li.matrixInOptimizer(mat2);
		}
		return mat2;
	}
	
	private boolean[][] notifyBreakFlagListener(boolean[][] brkFlag){
		for (MatrixListener li : matrixListeners) {
			brkFlag = li.initializeBreakFlag(brkFlag);
		}
		return brkFlag;
	}

	public void addMatrixListener(MatrixListener li){
		matrixListeners.add(li);
	}
	
	
	/**
	 * On input, mat[i][j] should give the score for aligning positions i and j.
	 * On output, mat[i][j] gives the maximum score possible for aligning 1..i 
	 * of protein 1 with 1..j of protein 2.
	 * 
	 * @param nSeq1 The length of protein 1 (mat.length)
	 * @param nSeq2 The length of protein 2 (mat[0].length)
	 * @param gapI gap initiation penalty
	 * @param gapE gap extension penalty
	 * @param isGlobal1 The alignment is global for protein 1
	 * @param isGlobal2 The alignment is global for protein 2
	 * @return The maximum score
	 */
	private double dpAlign(int nSeq1, int nSeq2, double gapI, double gapE, 
			boolean isGlobal1, boolean isGlobal2) {

		// isGlobal1,isGlobal2 are always false...

		int i, j, is, js, iMax, jMax, k;
		boolean ge=(gapE!=0.0?true:false);
		double sum, sum_ret, sum_brk;

		boolean[][] brk_flg=new boolean [nSeq1][nSeq2];
		for(i=0; i<nSeq1; i++) brk_flg[i]=new boolean [nSeq2];

		brk_flg = notifyBreakFlagListener(brk_flg);
		
		// ge = true here...
		/*  
		  for(i=0; i<nSeq1; i++)
		   {
		     printf("\n");
		     for(j=0; j<nSeq2; j++)
		       {
			 printf("%4d", (int)(*(mat[i]+j)*10));
		       }
		   }
		 printf("\n\n\n");
		 */
		if(!ge)
		{
			for(i=nSeq1-1; i>=0; i--)
				for(j=nSeq2-1; j>=0; j--)
				{
					brk_flg[i][j]=false;
					if(j<nSeq2-1 && i<nSeq1-1) 
					{
						sum=mat[i+1][j+1];
					}
					else
					{
						sum=0.0;
						if((isGlobal1 && i!=nSeq1-1) || (isGlobal2 && j!=nSeq2-1)) 
							sum=-gapI;
					}
					if(j+1<nSeq2)
						for(k=i+2; k<nSeq1; k++)
						{
							if(mat[k][j+1]-gapI>sum)
								sum=mat[k][j+1]-gapI;
						}
					if(i+1<nSeq1)
						for(k=j+2; k<nSeq2; k++)
						{
							if(mat[i+1][k]-gapI>sum)
								sum=mat[i+1][k]-gapI;
						}
					sum+=mat[i][j];
					sum_brk=(isGlobal1?-gapI:0.0)+(isGlobal2?-gapI:0.0);
					if(sum<sum_brk) 
					{
						sum=sum_brk;
						brk_flg[i][j]=true;
						//System.out.println("break at: " + i + " " + j);
					}
					mat[i][j]=sum;
				}
		}
		else
		{
			for(i=nSeq1-1; i>=0; i--)
				for(j=nSeq2-1; j>=0; j--)
				{
					brk_flg[i][j]=false;
					if(j<nSeq2-1 && i<nSeq1-1) 
					{
						sum=mat[i+1][j+1];
					}
					else 
					{
						sum=0.0;
						if(isGlobal1 && i!=nSeq1-1) sum=-gapI-gapE*(nSeq1-i-1);
						if(isGlobal2 && j!=nSeq2-1) sum=-gapI-gapE*(nSeq2-j-1);
					}
					if(j+1<nSeq2)
						for(k=i+2; k<nSeq1; k++)
							if(mat[k][j+1]-gapI-gapE*(k-i-1)>sum)
								sum=mat[k][j+1]-gapI-gapE*(k-i-1);
					if(i+1<nSeq1)
						for(k=j+2; k<nSeq2; k++)
							if(mat[i+1][k]-gapI-gapE*(k-j-1)>sum)
								sum=mat[i+1][k]-gapI-gapE*(k-j-1);
					sum+=mat[i][j];
					sum_brk=(isGlobal1?(-gapI-gapE*(nSeq1-1-i)):0.0)+(isGlobal2?(-gapI-gapE*(nSeq2-1-j)):0.0);
					if(sum<sum_brk) 
					{
						sum=sum_brk;
						brk_flg[i][j]=true;						 
					}
					mat[i][j]=sum;
				}
		}

		//		if (debug ){
		//			ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		//			JFrame frame = new JFrame("CE alignment matrix in dpAlign " );
		//			frame.addWindowListener(new WindowAdapter(){
		//				public void windowClosing(WindowEvent e){
		//					JFrame f = (JFrame) e.getSource();
		//					f.setVisible(false);
		//					f.dispose();
		//				}
		//
		//			});
		//
		//			smp.getMatrixPanel().setScalevalue(100);
		//			Matrix mx = (Matrix) new Matrix(mat).clone();
		//			smp.setMatrix(mx);
		//
		//			frame.getContentPane().add(smp);
		//
		//			frame.pack();
		//			frame.setVisible(true);
		//		}	

		/*  
		 for(i=0; i<nSeq1; i++)
		   {
		     printf("\n");
		     for(j=0; j<nSeq2; j++)
		       {
			 printf("%4d", (int)(*(mat[i]+j)*10));
		       }
		   }
		 printf("\n\n\n");
		 for(i=0; i<nSeq1; i++)
		   {
		     printf("\n");
		     for(j=0; j<nSeq2; j++)
		       {
			 printf("%4d", (int)(*(brk_flg[i]+j)));
		       }
		   }
		 // exit(0);
		 */

		is=0; js=0; lcmp=0;
		// no nc-end penalty - begin 
		sum_ret=mat[0][0];

		// look for the highest score in mat[i][j]
		for(i=0; i<nSeq1; i++)
			for(j=0; j<nSeq2; j++)
			{
				if(i==0 && j==0) continue;
				sum=mat[i][j];
				if(isGlobal1) sum+=-gapI-gapE*i;
				if(isGlobal2) sum+=-gapI-gapE*j;
				if(sum>sum_ret) 
				{
					sum_ret=sum;
					is=i; js=j;
				}
			}

		//System.out.println("start at " + is + "  " + js);
		//for(k=0; k<is; k++) align1[k]=-1;
		//for(k=0; k<js; k++) align2[k]=-1;
		// no nc-end penalty - end 

		for(i=is, j=js; i<nSeq1 && j<nSeq2; i++, j++)
		{
			iMax=i; jMax=j;
			sum=mat[i][j];
			if(!ge)
			{
				for(k=i+1; k<nSeq1; k++)
					if(mat[k][j]-gapI>sum)
					{
						iMax=k; jMax=j;
						sum=mat[k][j]-gapI;
					}

				for(k=j+1; k<nSeq2; k++)
					if(mat[i][k]-gapI>sum)
					{
						iMax=i; jMax=k;
						sum=mat[i][k]-gapI;
					}
			}
			else
			{
				for(k=i+1; k<nSeq1; k++)
					if(mat[k][j]-gapI-gapE*(k-i)>sum)
					{
						//System.out.println("gap1 " + k + " " + j + " " + sum + "<" +(mat[k][j]-gapI-gapE*(k-i)));
						iMax=k; jMax=j;
						sum=mat[k][j]-gapI-gapE*(k-i);
					}

				for(k=j+1; k<nSeq2; k++)
					if(mat[i][k]-gapI-gapE*(k-j)>sum)
					{
						//System.out.println("gap2 " + i + " " + k + " " + sum + "<"+ (mat[i][k]-gapI-gapE*(k-j)));
						iMax=i; jMax=k;
						sum=mat[i][k]-gapI-gapE*(k-j);
					}
			}

			//if ( i != iMax || j != jMax )
			//	System.out.println("FOUND GAP AT: " + i+ " " + iMax + " " + j + " " + jMax);

			//System.out.println(" iMax " + iMax + " jMax " +  jMax);
			// set the gap positions:
			//printf("%d %d\n", iMax, jMax);


			for(k=i; k<iMax; k++, i++) {
				align_se1[lcmp]=k;
				align_se2[lcmp]=-1;


				lcmp++;
			}

			for(k=j; k<jMax; k++, j++) {
				align_se1[lcmp]=-1;
				align_se2[lcmp]=k;


				lcmp++;
			}


			align_se1[lcmp]=iMax;
			align_se2[lcmp]=jMax; 
			lcmp++;

			if(brk_flg[i][j]) {
				//System.out.println("hit break flag at: " + i + "  " + j + " sum " + sum_ret + " lcmp " + lcmp);				
				break;

			}
		}


		return sum_ret;
	}




	private void rot_mol(Atom[] caA, Atom[] caB, int nse2, Matrix m , Atom shift) throws StructureException{



		for(int l=0; l<nse2; l++) {
			Atom a = caA[l];
			Group g = (Group)a.getGroup().clone();
			//Group g = (Group)a.getParent();

			Calc.rotate( g, m);
			Calc.shift(  g, shift);
			caB[l] = g.getAtom(a.getName());
		}
	}

	//void rot_mol(XYZ *molA, XYZ *molB, int nAtom, double d_[20] ) {
	//			  double dx, dy, dz;
	//			  for(int l=0; l<nAtom; l++) {
	//			    if(molA[l].X<1e10) {
	//			      dx=molA[l].X; 
	//			      dy=molA[l].Y;
	//			      dz=molA[l].Z;
	//			      molB[l].X=dx*d_[0]+dy*d_[1]+dz*d_[2]+d_[9];
	//			      molB[l].Y=dx*d_[3]+dy*d_[4]+dz*d_[5]+d_[10];
	//			      molB[l].Z=dx*d_[6]+dy*d_[7]+dz*d_[8]+d_[11];
	//			    }  
	//			    else {
	//			      molB[l]=molA[l];
	//			    }
	//			  }
	//			}
	//		


	/** superimpose and get rmsd
	 * 
	 * @param pro1
	 * @param pro2
	 * @param strLen
	 * @param storeTransform
	 * @param show Ignored. Formerly displayed the superposition with jmol.
	 * @return RMSD
	 * @throws StructureException
	 * @deprecated Use {@link #calc_rmsd(Atom[],Atom[],int,boolean)} instead
	 */
	@Deprecated
	public double calc_rmsd(Atom[] pro1, Atom[] pro2, int strLen, boolean storeTransform, boolean show) throws StructureException {
		return calc_rmsd(pro1, pro2, strLen, storeTransform);
	}

	/** superimpose and get rmsd
	 * 
	 * @param pro1
	 * @param pro2
	 * @param strLen Number of atoms from pro1 and pro2 to use
	 * @param storeTransform Store rotation and shift matrices locally
	 * @return RMSD
	 * @throws StructureException
	 */
	public double calc_rmsd(Atom[] pro1, Atom[] pro2, int strLen, boolean storeTransform) throws StructureException {

		Atom[] cod1 = getAtoms(pro1,  strLen,false);
		Atom[] cod2 = getAtoms(pro2,  strLen,true);

		assert(cod1.length == cod2.length);
		SVDSuperimposer svd = new SVDSuperimposer(cod1, cod2);

		Matrix matrix = svd.getRotation();
		Atom shift = svd.getTranslation();

		if ( storeTransform) {
			r=matrix;
			t = shift;
		}
		for (Atom a : cod2){
			Calc.rotate(a.getGroup(), matrix);
			Calc.shift(a.getGroup(),  shift);
		}
		return SVDSuperimposer.getRMS(cod1, cod2);

	}

	/**
	 * Copies the first length atoms from the input array
	 * @param ca The array to copy
	 * @param length the number of atoms to copy
	 * @param clone If true, preform a deep copy, cloning the underlying Groups
	 * @return An array with the first length items of ca, possibly cloning the Atoms.
	 * @throws StructureException
	 */
	private Atom[] getAtoms(Atom[] ca,  int length, boolean clone) throws StructureException{

		List<Atom> atoms = new ArrayList<Atom>();
		for ( int i = 0 ; i < length ; i++){

			Atom a;
			if ( clone ){
				Group g = (Group)ca[i].getGroup().clone();
				a = g.getAtom(ca[i].getName());
			}
			else {
				a = ca[i];
			}
			atoms.add(a);
		}
		return atoms.toArray(new Atom[atoms.size()]);
	}



	private void noBestTrace(){

		if(isPrint) {
			timeEnd = System.currentTimeMillis();
			long time_q=(timeEnd-timeStart);

			String msg = String.format("size=0 time=%d comb=%d\n", (int)(time_q), nTraces);
			System.out.println(msg);
		}
	}



	private double zToP(double z) {
		int ind=(int)(z/0.1);
		if(ind<0) ind=0;
		if(ind>149) ind=149;
		return(tableZtoP[ind]);
	}
	///////////////////////////////////////////////////////////////////////////
	private		double pToZ(double p) {
		int ind=(int)(-Math.log10(p)*3.0);
		if(ind<0) ind=0;
		if(ind>149) ind=149;
		return(tablePtoZ[ind]);
	}
	///////////////////////////////////////////////////////////////////////////
	private	double zByZ(double z1, double z2) {
		double p1=zToP(z1);
		double p2=zToP(z2);
		return(pToZ(p1*p2));
	}

	protected double zStrAlign(int winSize, int nTrace, double score, int nGaps) {
		double z1=zScore(winSize, nTrace, score);
		double z2=zGaps(winSize, nTrace, nGaps);
		return(zByZ(z1, z2));
	}

	double zScore(int winSize, int nTrace, double score) {

		if(winSize==8) {

			if(nTrace<1) return(0.0);

			double scoreAv_, scoreSd_;
			if(nTrace<21) {
				scoreAv_=scoreAv8[nTrace-1];
				scoreSd_=scoreSd8[nTrace-1];
			}
			else {
				scoreAv_=0.209874*nTrace+2.944714;
				scoreSd_=0.039487*nTrace+0.675735;
			}
			if(score>scoreAv_) return(0.0);
			return((scoreAv_-score)/scoreSd_);
		}

		if(winSize==6) {

			if(nTrace<1) return(0.0);

			double scoreAv_, scoreSd_;
			if(nTrace<21) {
				scoreAv_=scoreAv6[nTrace-1];
				scoreSd_=scoreSd6[nTrace-1];
			}
			else {
				scoreAv_=0.198534*nTrace+2.636477;
				scoreSd_=0.040922*nTrace+0.715636;
			}
			if(score>scoreAv_) return(0.0);
			return((scoreAv_-score)/scoreSd_);
		}

		return(0.0);

	}
	///////////////////////////////////////////////////////////////////////////
	double zGaps(int winSize, int nTrace, int nGaps) {

		if(nTrace<2) return(0.0);
		double scoreAv_, scoreSd_;

		if(winSize==8) {
			if(nTrace<21) {
				scoreAv_=gapsAv8[nTrace-1];
				scoreSd_=gapsSd8[nTrace-1];
			}
			else {
				scoreAv_=14.949173*nTrace-14.581193;
				scoreSd_=2.045067*nTrace+13.191095;
			}
			if(nGaps>scoreAv_) return(0.0);
			return((scoreAv_-nGaps)/scoreSd_);
		}

		if(winSize==6) {
			if(nTrace<21) {
				scoreAv_=gapsAv6[nTrace-1];
				scoreSd_=gapsSd6[nTrace-1];
			}
			else {
				scoreAv_=13.574490*nTrace-13.977223;
				scoreSd_=1.719977*nTrace+19.615014;
			}
			if(nGaps>scoreAv_) return(0.0);
			return((scoreAv_-nGaps)/scoreSd_);
		}

		return(0.0);
	}

	private static final double scoreAv8[]={2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74,
		4.99, 5.22, 5.46, 5.70, 5.94, 6.13, 6.36, 6.52, 6.68, 6.91};
	private static final double scoreSd8[]={1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04,
		1.08, 1.10, 1.15, 1.19, 1.23, 1.25, 1.32, 1.34, 1.36, 1.45};
	private static final double gapsAv8[]={0.00, 11.50, 23.32, 35.95, 49.02, 62.44, 76.28, 90.26, 
		104.86, 119.97, 134.86, 150.54, 164.86, 179.57, 194.39, 
		209.38, 224.74, 238.96, 253.72, 270.79};
	private static final double gapsSd8[]={0.00, 9.88, 14.34, 17.99, 21.10, 23.89, 26.55, 29.00, 31.11,
		33.10, 35.02, 36.03, 37.19, 38.82, 41.04, 43.35, 45.45, 
		48.41, 50.87, 52.27};
	private static final double scoreAv6[]={1.98, 1.97, 2.22, 2.54, 2.87, 3.18, 3.48, 3.77, 4.05, 4.31,
		4.57, 4.82, 5.03, 5.24, 5.43, 5.64, 5.82, 6.02, 6.21, 6.42};
	private static final double scoreSd6[]={1.15, 0.73, 0.63, 0.64, 0.71, 0.80, 0.87, 0.95, 1.01, 1.07,
		1.13, 1.19, 1.22, 1.25, 1.28, 1.32, 1.35, 1.39, 1.45, 1.50};
	private static final double gapsAv6[]={0.00, 10.12, 20.25, 31.29, 42.95, 55.20, 67.53, 80.15, 
		93.30, 106.47, 120.52, 134.38, 148.59, 162.58, 176.64, 
		191.23, 204.12, 218.64, 231.82, 243.43};
	private static final double gapsSd6[]={0.00, 9.80, 14.44, 18.14, 21.35, 24.37, 27.00, 29.68, 32.22,
		34.37, 36.65, 38.63, 40.31, 42.16, 43.78, 44.98, 47.08, 
		49.09, 50.78, 52.15};


	private static final double tableZtoP[]={
		1.0, 9.20e-01,8.41e-01,7.64e-01,6.89e-01,6.17e-01,5.49e-01,4.84e-01,4.24e-01,3.68e-01,
		3.17e-01,2.71e-01,2.30e-01,1.94e-01,1.62e-01,1.34e-01,1.10e-01,8.91e-02,7.19e-02,5.74e-02,
		4.55e-02,3.57e-02,2.78e-02,2.14e-02,1.64e-02,1.24e-02,9.32e-03,6.93e-03,5.11e-03,3.73e-03,
		2.70e-03,1.94e-03,1.37e-03,9.67e-04,6.74e-04,4.65e-04,3.18e-04,2.16e-04,1.45e-04,9.62e-05,
		6.33e-05,4.13e-05,2.67e-05,1.71e-05,1.08e-05,6.80e-06,4.22e-06,2.60e-06,1.59e-06,9.58e-07,
		5.73e-07,3.40e-07,1.99e-07,1.16e-07,6.66e-08,3.80e-08,2.14e-08,1.20e-08,6.63e-09,3.64e-09,
		1.97e-09,1.06e-09,5.65e-10,2.98e-10,1.55e-10,8.03e-11,4.11e-11,2.08e-11,1.05e-11,5.20e-12,
		2.56e-12,1.25e-12,6.02e-13,2.88e-13,1.36e-13,6.38e-14,2.96e-14,1.36e-14,6.19e-15,2.79e-15,
		1.24e-15,5.50e-16,2.40e-16,1.04e-16,4.46e-17,1.90e-17,7.97e-18,3.32e-18,1.37e-18,5.58e-19,
		2.26e-19,9.03e-20,3.58e-20,1.40e-20,5.46e-21,2.10e-21,7.99e-22,3.02e-22,1.13e-22,4.16e-23,
		1.52e-23,5.52e-24,1.98e-24,7.05e-25,2.48e-25,8.64e-26,2.98e-26,1.02e-26,3.44e-27,1.15e-27,
		3.82e-28,1.25e-28,4.08e-29,1.31e-29,4.18e-30,1.32e-30,4.12e-31,1.27e-31,3.90e-32,1.18e-32,
		3.55e-33,1.06e-33,3.11e-34,9.06e-35,2.61e-35,7.47e-36,2.11e-36,5.91e-37,1.64e-37,4.50e-38,
		1.22e-38,3.29e-39,8.77e-40,2.31e-40,6.05e-41,1.56e-41,4.00e-42,1.02e-42,2.55e-43,6.33e-44,
		1.56e-44,3.80e-45,9.16e-46,2.19e-46,5.17e-47,1.21e-47,2.81e-48,6.45e-49,1.46e-49,3.30e-50};
	private static final double tablePtoZ[]={
		0.00,0.73,1.24,1.64,1.99,2.30,2.58,2.83,3.07,3.29,
		3.50,3.70,3.89,4.07,4.25,4.42,4.58,4.74,4.89,5.04,
		5.19,5.33,5.46,5.60,5.73,5.86,5.99,6.11,6.23,6.35,
		6.47,6.58,6.70,6.81,6.92,7.02,7.13,7.24,7.34,7.44,
		7.54,7.64,7.74,7.84,7.93,8.03,8.12,8.21,8.30,8.40,
		8.49,8.57,8.66,8.75,8.84,8.92,9.01,9.09,9.17,9.25,
		9.34,9.42,9.50,9.58,9.66,9.73,9.81,9.89,9.97,10.04,
		10.12,10.19,10.27,10.34,10.41,10.49,10.56,10.63,10.70,10.77,
		10.84,10.91,10.98,11.05,11.12,11.19,11.26,11.32,11.39,11.46,
		11.52,11.59,11.66,11.72,11.79,11.85,11.91,11.98,12.04,12.10,
		12.17,12.23,12.29,12.35,12.42,12.48,12.54,12.60,12.66,12.72,
		12.78,12.84,12.90,12.96,13.02,13.07,13.13,13.19,13.25,13.31,
		13.36,13.42,13.48,13.53,13.59,13.65,13.70,13.76,13.81,13.87,
		13.92,13.98,14.03,14.09,14.14,14.19,14.25,14.30,14.35,14.41,
		14.46,14.51,14.57,14.62,14.67,14.72,14.77,14.83,14.88,14.93};

	/** copy data from this class into AFPChain container object.
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 */
	 public void convertAfpChain(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {

		 afpChain.setBlockNum(1);
		 //afpChain.setAlignScore(z);
		 Matrix[] m ;

		 if ( r != null ) {
			 m = new Matrix[1];
			 m[0] = r;
		 } else  {
			 m = new Matrix[0];
		 }

		 Atom[] as ;
		 if ( t != null) {
			 as = new Atom[1];
			 as[0] = t;
		 } else {
			 as = new Atom[0];
		 }

		 afpChain.setBlockRotationMatrix(m);
		 afpChain.setBlockShiftVector(as);

		 int nse1 = ca1.length;
		 int nse2 = ca2.length;
		 //System.out.println("dist1 :" + dist1.length + " " + dist2.length);

		 if ( nse1 > 0 && dist1.length > 0 )
			 afpChain.setDisTable1(new Matrix(dist1));
		 else 
			 afpChain.setDisTable1 (Matrix.identity(3, 3));
		 if ( nse2 > 0 && dist2.length > 0 )
			 afpChain.setDisTable2(new Matrix(dist2));
		 else
			 afpChain.setDisTable2(Matrix.identity(3, 3));


		 char[] alnseq1 = new char[nse1+nse2+1];
		 char[] alnseq2 = new char[nse1+nse2+1] ;
		 char[] alnsymb = new char[nse1+nse2+1];

		 int[][][] optAln = new int[1][2][nAtom];
		 afpChain.setOptAln(optAln);

		 int pos = 0;
		 int nrIdent = 0;
		 int nrSim   = 0;
		 for(int ia=0; ia<lcmp; ia++) {

			 // no gap
			 if(align_se1[ia]!=-1 && align_se2[ia]!=-1) {
				 //System.out.println("ia " + ia + " pos " + pos + " "  + align_se1[ia] + " " + align_se2[ia]);
				 optAln[0][0][pos] = align_se1[ia];
				 optAln[0][1][pos] = align_se2[ia];

				 char l1 = getOneLetter(ca1[align_se1[ia]].getGroup());
				 char l2 = getOneLetter(ca2[align_se2[ia]].getGroup());

				 alnseq1[ia] = Character.toUpperCase(l1);
				 alnseq2[ia] = Character.toUpperCase(l2);
				 alnsymb[ia] = ' ';
				 if ( l1 == l2) {					
					 nrIdent++;
					 nrSim++;
					 alnsymb[ia] = '|';
				 } else if ( AFPAlignmentDisplay.aaScore(l1, l2) > 0){
					 nrSim++;
					 alnsymb[ia] = ':';
				 }

				 pos++;

			 } else {
				 // there is a gap at this position
				 alnsymb[ia] = ' ';
				 if (align_se1[ia] == -1 ) {
					 alnseq1[ia] = '-';
				 } else {
					 char l1 = getOneLetter(ca1[align_se1[ia]].getGroup());
					 alnseq1[ia] = Character.toUpperCase(l1);
				 }
				 if ( align_se2[ia] == -1 ) {
					 alnseq2[ia] = '-';
				 } else {
					 char l2 = getOneLetter(ca2[align_se2[ia]].getGroup());
					 alnseq2[ia] = Character.toUpperCase(l2);
				 }

			 }
		 }


		 afpChain.setAlnseq1(alnseq1);
		 afpChain.setAlnseq2(alnseq2);
		 afpChain.setAlnsymb(alnsymb);


		 // CE uses the aligned pairs as reference not the whole alignment including gaps...
		 if ( pos > 0) {
			 afpChain.setIdentity(nrIdent*1.0/pos);
			 afpChain.setSimilarity(nrSim*1.0/pos);
		 } else {
			 afpChain.setIdentity(0);
			 afpChain.setSimilarity(0);
		 }

		 //AFPAlignmentDisplay.getAlign( afpChain,ca1,ca2);

	 }

	 public int getnAtom() {
		 return nAtom;
	 }



	 public void setnAtom(int nAtom) {
		 this.nAtom = nAtom;
	 }



	 public int getLcmp() {
		 return lcmp;
	 }



	 public void setLcmp(int lcmp) {
		 this.lcmp = lcmp;
	 }



	 public int[] getAlign_se1() {
		 return align_se1;
	 }



	 public void setAlign_se1(int[] alignSe1) {
		 align_se1 = alignSe1;
	 }



	 public int[] getAlign_se2() {
		 return align_se2;
	 }



	 public void setAlign_se2(int[] alignSe2) {
		 align_se2 = alignSe2;
	 }

	 /**
	  * Caution: this matrix is overwriten with very different data at several
	  * points in the alignment algorithm. After 
	  * {@link #initSumOfDistances(int, int, int, int, Atom[], Atom[]) initSumOfDistances}
	  * is run, this will hold the distance matrix between AFPs.
	  * @return mat
	  */
	 public double[][] getMatMatrix() {
		 return mat;
	 }
	 
	 public void setMatMatrix(double[][] matrix){
		 mat = matrix;
	 }

	 /**
	  * Gets the rotation matrix from the last call to 
	  * {@link #calc_rmsd(Atom[], Atom[], int, boolean) calc_rmsd}.
	  * @return The rotatiokn matrix
	  */
	 public Matrix getRotationMatrix() {
		 return r;
	 }

	 /**
	  * Gets the shift from the last call to 
	  * {@link #calc_rmsd(Atom[], Atom[], int, boolean) calc_rmsd}.
	  * @return The shift
	  */
	 public Atom getShift() {
		 return t;
	 }

	public double[][] getDist1() {
		return dist1;
	}

	public void setDist1(double[][] dist1) {
		this.dist1 = dist1;
	}

	public double[][] getDist2() {
		return dist2;
	}

	public void setDist2(double[][] dist2) {
		this.dist2 = dist2;
	}
	 
	 
}
