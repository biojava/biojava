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
package org.biojava.nbio.tornadovm.bench;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.GuideTree;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.tornadovm.TornadoSupport;
import org.biojava.nbio.tornadovm.align.TornadoAlignments;

/**
 * All-pairs alignment scores: BioJava (thread pool) vs TornadoVM.
 * Usage: SwBench fasta:FILE | random:N:MINLEN:MAXLEN ...; property bench.types (default GLOBAL,LOCAL).
 */
public class SwBench {

	public static void main(String[] args) throws Exception {
		SubstitutionMatrix<AminoAcidCompound> blosum = SubstitutionMatrixHelper.getBlosum62();
		SimpleGapPenalty gaps = new SimpleGapPenalty();
		String[] types = System.getProperty("bench.types", "GLOBAL,LOCAL").split(",");
		boolean skipCpu = Boolean.getBoolean("bench.skipcpu");
		System.out.printf("threads in BioJava pool: %d%n", ConcurrencyTools.getThreadPool().getMaximumPoolSize());
		System.out.printf("%-26s %-6s %9s %8s | %11s %11s %10s | %9s %9s %7s | %s%n", "set", "type", "pairs", "Gcells",
				"biojava ms", "leanCPU ms", "gpu ms", "vsBioJava", "vsLeanCPU", "GCUPS", "parity");
		for (String spec : args) {
			List<ProteinSequence> seqs = load(spec);
			long cells = 0;
			for (int i = 0; i < seqs.size(); i++)
				for (int j = i + 1; j < seqs.size(); j++)
					cells += (long) seqs.get(i).getLength() * seqs.get(j).getLength();
			String dump = System.getProperty("bench.dump");
			String tag = spec.replaceAll("[^A-Za-z0-9]+", "_");
			if (dump != null) {
				try (java.io.PrintWriter w = new java.io.PrintWriter(new File(dump, tag + ".fasta"))) {
					for (int i = 0; i < seqs.size(); i++) {
						w.println(">s" + i);
						w.println(seqs.get(i).getSequenceAsString());
					}
				}
			}
			for (String t : types) {
				PairwiseSequenceScorerType type = PairwiseSequenceScorerType.valueOf(t);
				double[] ref = null;
				double cpu = Double.NaN;
				if (!skipCpu) {
					long t0 = System.nanoTime();
					ref = Alignments.getAllPairsScores(seqs, type, gaps, blosum);
					cpu = (System.nanoTime() - t0) / 1e6;
				}
				int[][] codes = new int[seqs.size()][];
				java.util.Map<AminoAcidCompound, Integer> idx = new java.util.HashMap<>();
				java.util.List<AminoAcidCompound> alpha = new ArrayList<>();
				for (int i = 0; i < seqs.size(); i++) {
					codes[i] = new int[seqs.get(i).getLength()];
					int pos = 0;
					for (AminoAcidCompound c : seqs.get(i)) {
						Integer v = idx.get(c);
						if (v == null) { v = alpha.size(); idx.put(c, v); alpha.add(c); }
						codes[i][pos++] = v;
					}
				}
				int kk = alpha.size();
				int[] subs = new int[kk * kk];
				for (int a = 0; a < kk; a++) for (int b = 0; b < kk; b++) subs[a * kk + b] = blosum.getValue(alpha.get(a), alpha.get(b));
				double lean = Double.MAX_VALUE;
				int[] leanScores = null;
				for (int r = 0; r < 3; r++) {
					long t0 = System.nanoTime();
					leanScores = org.biojava.nbio.tornadovm.align.CpuAlignmentScores.allPairs(codes, subs, kk, gaps.getOpenPenalty(), gaps.getExtensionPenalty(), type == PairwiseSequenceScorerType.LOCAL);
					lean = Math.min(lean, (System.nanoTime() - t0) / 1e6);
				}
				System.setProperty(TornadoSupport.PROPERTY, "force");
				double gpu = Double.MAX_VALUE;
				double[] got = null;
				for (int r = 0; r < 3; r++) {
					long t0 = System.nanoTime();
					got = TornadoAlignments.getAllPairsScores(seqs, type, gaps, blosum);
					gpu = Math.min(gpu, (System.nanoTime() - t0) / 1e6);
				}
				System.clearProperty(TornadoSupport.PROPERTY);
				if (dump != null) {
					try (java.io.PrintWriter w = new java.io.PrintWriter(new File(dump, tag + "." + t + ".scores"))) {
						for (double v : got) w.println((int) v);
					}
				}
				String parity = "biojava skipped";
				int leanDiff = 0;
				for (int k = 0; k < got.length; k++) if (leanScores[k] != got[k]) leanDiff++;
				if (leanDiff > 0) parity = "LEAN CPU DIFFERS on " + leanDiff;
				if (ref != null) {
					int diff = 0;
					for (int k = 0; k < ref.length; k++)
						if (ref[k] != got[k]) diff++;
					parity = diff == 0 ? "IDENTICAL (" + ref.length + " scores)" : diff + " of " + ref.length + " DIFFER";
				}
				System.out.printf("%-26s %-6s %9d %8.2f | %11.1f %11.1f %10.1f | %8.1fx %8.1fx %7.1f | %s%n", spec, t,
						got.length, cells / 1e9, cpu, lean, gpu, cpu / gpu, lean / gpu, cells / gpu / 1e6, parity);
			}
			if (Boolean.getBoolean("bench.tree")) {
				List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> cpuScorers = Alignments
						.getAllPairsScorers(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, blosum);
				Alignments.runPairwiseScorers(cpuScorers);
				System.setProperty(TornadoSupport.PROPERTY, "force");
				List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> gpuScorers = TornadoAlignments
						.getAllPairsScorers(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, blosum);
				System.clearProperty(TornadoSupport.PROPERTY);
				String a = new GuideTree<>(seqs, cpuScorers).toString();
				String b = new GuideTree<>(seqs, gpuScorers).toString();
				System.out.printf("  guide tree (GLOBAL scores) for %s: %s%n", spec, a.equals(b) ? "IDENTICAL newick" : "DIFFERENT");
			}
		}
		ConcurrencyTools.shutdown();
	}

	static List<ProteinSequence> load(String spec) throws Exception {
		String[] p = spec.split(":");
		if (p[0].equals("fasta")) {
			// ungapped sequences (the file may be an alignment)
			List<ProteinSequence> seqs = new ArrayList<>();
			for (ProteinSequence s : FastaReaderHelper.readFastaProteinSequence(new File(p[1])).values()) {
				seqs.add(new ProteinSequence(s.getSequenceAsString().replace("-", "")));
			}
			return seqs;
		}
		int n = Integer.parseInt(p[1]), minLen = Integer.parseInt(p[2]), maxLen = Integer.parseInt(p[3]);
		String aa = "ACDEFGHIKLMNPQRSTVWY";
		Random rnd = new Random(42);
		List<ProteinSequence> seqs = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			int len = minLen + rnd.nextInt(maxLen - minLen + 1);
			StringBuilder sb = new StringBuilder(len);
			for (int k = 0; k < len; k++) sb.append(aa.charAt(rnd.nextInt(aa.length())));
			seqs.add(new ProteinSequence(sb.toString()));
		}
		return seqs;
	}
}
