package org.biojava.bio.structure.align.benchmark.metrics;


import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;


public abstract class RMSDMetric {
	public static class Reference extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			try {

				List<Atom[]> structures = new ArrayList<Atom[]>(2);
				structures.add(ca1);
				structures.add(ca2);
				int[][] optAln = reference.getAlignmentMatrix(structures);

				// Create new arrays for the subset of atoms in the alignment.
				Atom[] ca1aligned = new Atom[reference.size()];
				Atom[] ca2aligned = new Atom[reference.size()];
				int pos=0;
				for(int i=0;i<optAln[0].length;i++) {
					ca1aligned[pos] = ca1[optAln[0][pos]];
					ca2aligned[pos] = (Atom) ca2[optAln[1][pos]].clone();
					pos++;

				}
				//Superimpose
				SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
				Matrix matrix = svd.getRotation();
				Atom shift = svd.getTranslation();

				for(Atom a : ca2aligned) {
					Calc.rotate(a, matrix);
					Calc.shift(a, shift);
				}

				return SVDSuperimposer.getRMS(ca1aligned, ca2aligned);
			} catch (StructureException e) {
				e.printStackTrace();
				return Double.NaN;
			}
		}

		@Override
		public String getName() {
			return "Ref_RMSD";
		}

		@Override
		public String format(double result) {
			return String.format("%.4f", result);
		}
	}

	public static class Alignment extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			// check if already set in the chain
			double optRMSD = align.getTotalRmsdOpt();
			if( optRMSD > 0 ) {
				return optRMSD;
			}
			// Create new arrays for the subset of atoms in the alignment.
			Atom[] ca1aligned = new Atom[align.getOptLength()];
			Atom[] ca2aligned = new Atom[align.getOptLength()];
			int pos=0;
			int[] blockLens = align.getOptLen();
			int[][][] optAln = align.getOptAln();
			assert(align.getBlockNum() <= optAln.length);
			for(int block=0;block< align.getBlockNum();block++) {
				assert(blockLens[block] <= optAln[block][0].length);
				for(int i=0;i<blockLens[block];i++) {
					ca1aligned[pos] = ca1[optAln[block][0][i]];
					ca2aligned[pos] = (Atom) ca2[optAln[block][1][i]].clone();
					pos++;
				}
			}

			try {

				//Superimpose
				SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
				Matrix matrix = svd.getRotation();
				Atom shift = svd.getTranslation();

				for(Atom a : ca2aligned) {
					Calc.rotate(a, matrix);
					Calc.shift(a, shift);
				}

				return SVDSuperimposer.getRMS(ca1aligned, ca2aligned);
			} catch (StructureException e) {
				e.printStackTrace();
				return Double.NaN;
			}
		}

		@Override
		public String getName() {
			return "Aln_RMSD";
		}

		@Override
		public String format(double result) {
			return String.format("%.4f", result);
		}
	}
	
	public static void main(String[] args) {
		String pdb1 = "121p";
		String pdb2 = "1htw.A";
		try {
			AtomCache cache = new AtomCache("/tmp/", true);
			
			Atom[] ca1 = cache.getAtoms(pdb1);
			Atom[] ca2 = cache.getAtoms(pdb2);

			CeMain ceMain;
			ceMain = (CeMain) StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			//((CeParameters)ceMain.getParameters()).setCheckCircular(true);
			((CeParameters)ceMain.getParameters()).setMaxGapSize(0);

			AFPChain align = ceMain.align(ca1, ca2);
			System.out.println(align.toCE(ca1, ca2));

			Atom[] ca1aligned = new Atom[align.getOptLength()];
			Atom[] ca2aligned = new Atom[align.getOptLength()];
			int pos=0;
			int[] blockLens = align.getOptLen();
			int[][][] optAln = align.getOptAln();
			for(int block=0;block< optAln.length;block++) {
				//for(int i=0;i<blockLens[block];i++) {
				assert(blockLens[block] <= optAln[block][0].length);
				for(int i=0;i<blockLens[block];i++) {
					ca1aligned[pos] = ca1[optAln[block][0][i]];
					ca2aligned[pos] = (Atom) ca2[optAln[block][1][i]].clone();
					pos++;
				}
			}

			//Superimpose
			SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
			Matrix matrix = svd.getRotation();
			Atom shift = svd.getTranslation();

			for(Atom a : ca2aligned) {
				Calc.rotate(a, matrix);
				Calc.shift(a, shift);
			}

			System.out.println( SVDSuperimposer.getRMS(ca1aligned, ca2aligned));

			System.out.println( matrix.toString());
			System.out.println( shift.toString());



			PrintWriter pw1 = new PrintWriter(new FileWriter("/Users/blivens/dev/bourne/benchmarks/last1.pdb"));
			Structure s1 = ca1[0].getParent().getParent().getParent();
			pw1.println(s1.toPDB());                   
			pw1.close();

			PrintWriter pw2 = new PrintWriter(new FileWriter("/Users/blivens/dev/bourne/benchmarks/last2.pdb"));
			Structure s2 = ca2[0].getParent().getParent().getParent();
			pw2.println(s2.toPDB());                   
			pw2.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}
}
