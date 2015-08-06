package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;

public class CESymmCalculator extends CECalculator {

	public static final int MIN_ANGLE = 20;

	Atom origin1 = null ;
	Atom origin2  = null;

	public CESymmCalculator(CeParameters params) {
		super(params);

	}

	@SuppressWarnings("unused")
	@Override
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

													double angle = checkAngle(mse1,mse2, ca2, ca2, winSize);
													if ( angle < MIN_ANGLE)
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

	/** do a SVN of the AFPs at positions mse1, mse2. Returns the rotation angle between the AFPs
	 * 
	 * @param mse1
	 * @param mse2
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	private double checkAngle(int mse1, int mse2, Atom[] ca1, Atom[] ca2, int winSize) {
		try {

			if ( origin1 == null)
				origin1 = Calc.getCentroid(ca1);
			if ( origin2 == null)
				origin2= Calc.getCentroid(ca2);

			int nse1 = ca1.length;
			int nse2 = ca2.length;


			int max1 = Math.min((nse1-mse1-1),winSize );
			int max2 = Math.min((nse2-mse2-1),winSize );
			
			int maxAtoms = Math.min(max1,max2);
			
			//System.out.println(max1 + "  " + max2 + " " + maxAtoms);
			Atom[] cod1 = new Atom[maxAtoms];
			Atom[] cod2 = new Atom[maxAtoms];
			
			assert(cod1.length == cod2.length);

			// this sums up over the distances of the fragments

			int is1 = mse1;
			for(int pos1=0; pos1< maxAtoms; pos1++){
				is1++;
				cod1[pos1] = (Atom)ca1[is1].clone();
				//System.out.println(pos1 + " " + is1 + " " + cod1[pos1]);
			}
			int is2 = mse2;
			for(int pos2=0; pos2<  maxAtoms; pos2++){
				is2++;
				cod2[pos2] = (Atom)ca2[is2].clone();
			}	

			SVDSuperimposer svd = new SVDSuperimposer(cod1, cod2);

			Matrix matrix = svd.getRotation();
			Atom shift = svd.getTranslation();

			for ( Atom a: cod2){
				Calc.rotate( a, matrix);
				Calc.shift(  a, shift);
			}
			
			Atom v1 = Calc.subtract(origin1, cod1[0]);
			Atom v2 = Calc.subtract(origin2, cod2[0]);

			double ang= Calc.angle(v1,v2);
			//System.out.println(v1 + " " + v2 + " " + ang);
			return ang;


		} catch (Exception e){
			e.printStackTrace();
			return 0d;
		}
	}

}
