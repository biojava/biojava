/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.benchmark.metrics.Metric;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;

/**
 * @author Spencer Bliven
 *
 */
public class AlignBenchmark {
	private AtomCache cache;
	private StructureAlignment aligner;
	private int maxLength; // Don't try alignments between proteins which are too long.
	private List<Metric> metrics;

	public AlignBenchmark(String cachePath, StructureAlignment method, int maxLength) {
		this.cache = new AtomCache(cachePath, true);
		this.aligner = method;
		this.maxLength = maxLength;
		this.metrics = AlignmentStats.defaultMetrics();
	}

	/* TODO write a test which makes sure that CeMain gives good performance against RIPC
	@Test
	public void runRIPCBenchmark() {
		String RIPCfile = "src/test/resources/align/benchmarks/RIPC.align";
		MultipleAlignmentParser parser = new RIPCParser(RIPCfile);
		runBenchmark(parser);
	}
	 */

	public void runBenchmark(Writer out, MultipleAlignmentParser parser) {

		try {
			AlignmentStatsIterator statsIt = new AlignmentStatsIterator(parser,
					aligner, metrics, this.maxLength );


			AlignmentStats.writeHeader(out, metrics);
			while(statsIt.hasNext()) {
				AlignmentStats stats = statsIt.next();
				stats.writeStats(out);
				out.flush();
			}

		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}

	/**
	 * 
	 * @author Spencer Bliven
	 *
	 */
	// TODO test that this works with MultipleAlignments of more than two proteins
	protected class AlignmentStatsIterator implements Iterator<AlignmentStats> {
		private List<Metric> metrics;

		private Iterator<MultipleAlignment> parser;
		private final int maxLength;
		private StructureAlignment aligner;

		private MultipleAlignment currentAlignment;
		private List<Atom[]> structures; //Atoms for each protein in the currentAlignment
		private int nextAlignLeft; //index to one protein in the currentAlignment
		private int nextAlignRight; //index to another protein in the currentAlignment

		public AlignmentStatsIterator( MultipleAlignmentParser parser,
				StructureAlignment aligner, List<Metric> metrics)
		throws StructureException, IOException
		{
			this(parser,aligner,metrics,-1);
		}

		public AlignmentStatsIterator( MultipleAlignmentParser parser,
				StructureAlignment aligner, List<Metric> metrics, int maxLength)
		throws StructureException, IOException
		{
			this.parser = parser.iterator();
			this.maxLength = maxLength;
			this.aligner = aligner;

			this.metrics = metrics;

			if(this.parser.hasNext()) {
				this.currentAlignment = this.parser.next();

				// Get alignment structures
				String[] pdbIDs = this.currentAlignment.getNames();
				this.structures= new ArrayList<Atom[]>(pdbIDs.length);

				for(int i=0;i<pdbIDs.length;i++) {
					String pdb1 = pdbIDs[i];

					Atom[] ca1 = cache.getAtoms(pdb1);

					if(ca1.length > this.maxLength && this.maxLength > 0) {
						System.err.format("Skipping large protein '%s' of length %d\n", pdb1, ca1.length);
						this.structures.add(null);
					} else {
						this.structures.add(ca1);
					}
				}

			} else {
				this.currentAlignment = null;
				this.structures = null;
			}
			this.nextAlignLeft=0;
			this.nextAlignRight=1;

		}

		/**
		 * @see java.util.Iterator#next()
		 */
		//@Override
		public AlignmentStats next() {
			while(true) { // Eventually returns or throws exception
				//String[] pdbIDs = this.currentAlignment.getNames();
				int alignLen = this.currentAlignment.getNames().length;

				while( nextAlignLeft<alignLen-1 ) {
					Atom[] ca1 = structures.get(nextAlignLeft);

					if(ca1 == null) {
						nextAlignLeft++;
						nextAlignRight = nextAlignLeft+1;
						continue;
					}

					while( nextAlignRight<alignLen) {
						Atom[] ca2 = structures.get(nextAlignRight);

						//Increment loop counter
						nextAlignRight++;

						if(ca2 == null) {
							continue;
						}

						//System.out.format("Comparing %s to %s\n",pdbIDs[i],pdbIDs[j]);
						AFPChain alignment;
						long elapsedTime;
						try {
							long startTime = System.currentTimeMillis();
							alignment = aligner.align(ca1,ca2);
							elapsedTime = System.currentTimeMillis() - startTime;
						} catch (StructureException e) {
							e.printStackTrace();
							return null;
						}
						//System.out.format("Done Comparing %s to %s\n",pdbIDs[i],pdbIDs[j]);

						Map<String,Object> metaData = new HashMap<String,Object>();
						metaData.put("alignmentTime",new Double(elapsedTime/1000.));
						return new AlignmentStats(this.metrics,this.currentAlignment, alignment, ca1, ca2, metaData);
					}

					nextAlignLeft++;
					nextAlignRight = nextAlignLeft+1;

				}


				// Nothing left for currentAlignment. Go on to next one.
				if(this.parser.hasNext()) {
					// Reset the currentAlignment
					this.currentAlignment = this.parser.next();
					nextAlignLeft = 0;
					nextAlignRight = 1;

					// Get alignment structures
					String[] pdbIDs = this.currentAlignment.getNames();

					this.structures= new ArrayList<Atom[]>(pdbIDs.length);

					for(int i=0;i<pdbIDs.length;i++) {
						String pdb1 = pdbIDs[i];

						Atom[] ca1;
						try {
							ca1 = cache.getAtoms(pdb1);
						} catch (StructureException e) {
							e.printStackTrace();
							return null;
						} catch (IOException e) {
							e.printStackTrace();
							return null;
						}

						if(ca1.length > this.maxLength && this.maxLength > 0) {
							System.err.format("Skipping large protein '%s' of length %d\n", pdb1, ca1.length);
							this.structures.add(null);
						} else {
							this.structures.add(ca1);
						}
					}
				}
				else {
					throw new IllegalStateException("No remaining items in iterator");
				}

			} // Return to loop start

		}

		/**
		 * @see java.util.Iterator#hasNext()
		 */
		//@Override
		public boolean hasNext() {
			int alignLen = this.currentAlignment.getNames().length;

			return this.parser.hasNext() ||
			nextAlignLeft >= alignLen-1 ||
			nextAlignRight<alignLen;
		}

		/**
		 * @see java.util.Iterator#remove()
		 */
		//@Override
		public void remove() {
			throw new UnsupportedOperationException("remove() not implemented");
		}

	}

	public static void main(String[] args) {
		//// Set parameters
		// TODO make these arguments
		int maxLength = 0; // Don't try alignments between proteins which are too long.
		String inFile = null;
		String outFile = null;
		MultipleAlignmentParser parser;
		StructureAlignment aligner;

		// Argument parsing
		String usage = "usage: AlignBenchmark parser alignment [infile [outfile [maxLength]]]\n" +
		"  parser:    Either RIPC or CPDB\n" +
		"  alignment: One of CE, CE0, CE-CP, CE-sidechains SmithWaterman, FATCAT-rigid, FATCAT-flexible\n" +
		"  infile:    If empty or omitted, gets alignments from biojava resource directory\n" +
		"  outfile:   If empty, omitted or \"-\", uses standard out\n" +
		"  maxLength: Longest protein to compare. 0 disables. Default 1000";

		if(args.length < 2) {
			System.err.println("Error: Too few args!");
			System.err.println(usage);
			System.exit(1);
			return;
		}

		// Do arg=0 (parser type) at the end
		int arg = 1;

		String alignerName = args[arg].toUpperCase();
		if(alignerName.substring(0, 2).equals("CE")) { // Different flavors of CE
			System.out.println("Using CE aligner");
			CeMain ceMain;
			try {
				ceMain = (CeMain) StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
				CeParameters param = (CeParameters)ceMain.getParameters();
				//CE0
				if(alignerName.equals("CE0")) {
					System.out.println("Setting Gap Size 0");
					param.setMaxGapSize(0);
				}
				//CE-CP
				else if(alignerName.matches("CE.?CP.*")) {
					System.out.println("Setting Gap Size 0");
					System.out.println("Checking for circular permutations");
					param.setMaxGapSize(0);
					param.setCheckCircular(true);
				}
				

				//-sidechains
				Pattern digitPat = Pattern.compile("(?:^|[^A-Z0-9])([0-9]+)(?:[^A-Z0-9]|$)");
				Matcher digitMatch = digitPat.matcher(alignerName);
				if( digitMatch.find() ) {
					System.out.println("Setting scoring mode to "+digitMatch.group(1));
					param.setScoringStrategy(Integer.parseInt(digitMatch.group(1)));
				}
				if(alignerName.contains("SIDECHAIN")) {
					System.out.println("Setting scoring mode to "+CeParameters.SIDE_CHAIN_SCORING);
					param.setScoringStrategy(CeParameters.SIDE_CHAIN_SCORING);
				}
				
			} catch (StructureException e) {
				e.printStackTrace();
				System.exit(1);
				return;
			}
			aligner = ceMain;
		}
		else if(alignerName.matches("SMITH-?WATERMAN|SW")) {
			System.out.println("Using Smith-Waterman aligner");

			try {

				aligner = StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName);
			} catch (StructureException e) {
				e.printStackTrace();
				System.exit(1);
				return;
			}
		}
		else if(alignerName.matches("J?FATCAT-?RIGID")) {
			System.out.println("Using FatCat-rigid aligner");

			try {
				aligner = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
			} catch (StructureException e) {
				e.printStackTrace();
				System.exit(1);
				return;
			}
		}
		else if(alignerName.matches("J?FATCAT-?FLEX.*")) {
			System.out.println("Using FatCat-flexible aligner");

			try {
				aligner = StructureAlignmentFactory.getAlgorithm(FatCatFlexible.algorithmName);
			} catch (StructureException e) {
				e.printStackTrace();
				System.exit(1);
				return;
			}
		}
		else {
			// Try matching against all known StructureAlignment algorithm names
			try {
				aligner = StructureAlignmentFactory.getAlgorithm(args[arg]);
				System.out.println("Using "+aligner.getAlgorithmName()+" aligner");
			} catch (StructureException e) {
				e.printStackTrace();
				System.exit(1);
				return;
			}
		}
		arg++;

		if(args.length > arg && !args[arg].equals("")) {
			inFile = args[arg];
		}
		arg++;

		if(args.length > arg && !args[arg].equals("")) {
			outFile = args[arg];
		}
		arg++;

		if(args.length > arg) {
			try {
				maxLength = Integer.parseInt(args[arg]);
			} catch(NumberFormatException e) {
				System.err.println("Unrecognized maximum length ("+args[arg]+").");
				System.err.println(usage);
				System.exit(1);
				return;
			}
		}
		arg++;

		if(args.length > arg) {
			System.err.println("Error: Too many args!");
			System.err.println(usage);
			System.exit(1);
			return;
		}
		arg++;

		String fileType = args[0];
		if(fileType.equalsIgnoreCase("RIPC")) {
			if(inFile == null )
				inFile = "src/test/resources/align/benchmarks/RIPC.align";
			//outFile = "/Users/blivens/dev/bourne/benchmarks/RIPC.stats";
			parser = new RIPCParser(inFile);
		}
		else if(fileType.equalsIgnoreCase("CPDB")) {
			if(inFile == null )
				inFile = "src/test/resources/align/benchmarks/CPDB_CPpairsAlignments_withCPSiteRefinement.txt";
			//outFile = "/Users/blivens/dev/bourne/benchmarks/CPDB.stats";
			parser = new CPDBParser(inFile);
		}
		else {
			System.err.println("Unrecognized parser. Choose from: RIPC, CPDB");
			System.err.println(usage);
			System.exit(1);
			return;
		}

		//// Done with parameters

		Writer out;
		out = new BufferedWriter(new OutputStreamWriter(System.out));
		if( outFile != null && !(outFile.equals("") || outFile.equals("-") )) {
			try {
				out = new BufferedWriter(new FileWriter(outFile));
			} catch (IOException e1) {
				e1.printStackTrace();
				System.err.println("Error writing file "+outFile);
				System.exit(1);
				return;
			}
		}

		AlignBenchmark bm = new AlignBenchmark("/tmp/",aligner,maxLength);
		bm.runBenchmark(out, parser);
	}

}
