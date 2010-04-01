/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.biojava.bio.structure.align.benchmark.MultipleAlignmentParser;

/**
 * Parses files formatted like the RIPC dataset (Repetition, Indels, Permutation,
 * Conformational variability), cited in:<br/>
 * <i>Mayr et al. Comparative analysis of protein structure alignments. BMC Struct Biol (2007) vol. 7 pp. 50</i>
 * <br/>and available online <a href="http://www.biomedcentral.com/content/supplementary/1472-6807-7-50-S5.txt">here</a>
 * or at src/test/resources/align/benchmarks/RIPC.align.
 * <p>
 * See {@link RIPCIterator} for info about the format.
 * @author Spencer Bliven
 *
 */
public class RIPCParser implements MultipleAlignmentParser
{
	private String filename;
	public RIPCParser(String filename) {
		this.filename=filename;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		String filename = "src/test/resources/align/benchmarks/RIPC.align";
		RIPCParser parser = new RIPCParser(filename);
		for(MultipleAlignment ma : parser) {
			System.out.println(ma.display());
		}
	}

	@Override
	public Iterator<MultipleAlignment> iterator()
	{
		try {
			return new RIPCIterator(filename);
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	

	private static final Pattern scopRegex = Pattern.compile("d(....)(.)(.)");
	/**
	 * 
	 * @param scopIDs
	 * @return
	 */
	public static String[] getPDBNames(String[] scopIDs) {
		String[] pdbNames = new String[scopIDs.length];
		for(int i=0;i<scopIDs.length;i++) {
			pdbNames[i] = getPDBName(scopIDs[i]);
		}
		return pdbNames;
	}
	/**
	 * Converts scop domain identifiers (eg 'd1lnlb1') into a PDB ID and chain
	 * (eg '1lnl.B').
	 * @param scopID
	 * @return the extracted pdbid and chain, or null if the scopID is malformed.
	 */
	public static String getPDBName(String scopID) {
		Matcher match = scopRegex.matcher(scopID);
		if(!match.matches()) {
			return null;
		}
		if(!match.group(2).equals("_")) {
			return match.group(1)+"."+match.group(2).toUpperCase();
		} else {
			return match.group(1);
		}
	}
	
	
	
	/**
	 * Parses a RIPC alignment file.
	 * <p>
	 * Format:<br/><pre>
	 * alignment := comment* labelPair resPair+ '\n'
	 * comment := '#' .* '\n' | '\n'
	 * labelPair := '#' (label1) '-' (label2)
	 * resPair := res ' ' res '\n'
	 * res := (3-letter amino acid code)'.'(pdb res number)'.'(insertion code)'.'(chain)</pre>
	 * Where outputs are listed in parentheses. Note that alignments must be separated by empty lines.
	 * <p>
	 * Example:<br/><pre>
	 * #d1hcy_2-d1lnlb1
	 * HIS.194._._	HIS.41._.B
	 * HIS.198._._	HIS.61._.B
	 * HIS.224._._	HIS.70._.B
	 * HIS.344._._	HIS.181._.B
	 *
	 * @author Spencer Bliven
	 *
	 */
	protected static class RIPCIterator implements Iterator<MultipleAlignment> {
		private BufferedReader ripc;
		private String[] nextLabels;

		private static final Pattern labelRegex =
			Pattern.compile("^\\#([^-]+)-([^-]+)$");
		private static final Pattern pairRegex = 
			Pattern.compile("^([A-Z]{3})\\.(\\d+)\\.(.)\\.(.)\\s+([A-Z]{3})\\.(\\d+)\\.(.)\\.(.)$");
		private static final Pattern commentRegex = 
			Pattern.compile("^(?:#.*|$)"); // labels are a subset of comments; check for labels first

		public RIPCIterator(String filename) throws IOException {
			ripc = new BufferedReader(new FileReader(filename));
			assert(ripc.markSupported());
			nextLabels = null;
			skipComments();
		}
		
		@Override
		public MultipleAlignment next()
		{
			List<List<PDBResidue>> residues = new ArrayList<List<PDBResidue>>(2);
			residues.add(new LinkedList<PDBResidue>());
			residues.add(new LinkedList<PDBResidue>());


			String line;
			try {
				line = ripc.readLine();

				while( line!=null ) {
					line = line.trim();
					Matcher pair = pairRegex.matcher(line);
					if(pair.matches()) {
						String aa1 = pair.group(1);
						String aa2 = pair.group(5);

						String pdb1 = pair.group(2);
						String pdb2 = pair.group(6);

						//Add insertion codes
						if( !pair.group(3).equals("_") ) {
							pdb1 += pair.group(3);
						}
						if( !pair.group(7).equals("_") ) {
							pdb1 += pair.group(7);
						}

						String chain1 = pair.group(4);
						String chain2 = pair.group(8);


						residues.get(0).add(new PDBResidue(pdb1,chain1,aa1));
						residues.get(1).add(new PDBResidue(pdb2,chain2,aa2));
					}
					else {
						Matcher labels=labelRegex.matcher(line);
						if(labels.matches()) {
							MultipleAlignment m = new MultipleAlignment(RIPCParser.getPDBNames(nextLabels),residues);

							nextLabels[0] = labels.group(1);
							nextLabels[1] = labels.group(2);

							return m;
						}
						else if(commentRegex.matcher(line).matches()) {
							// ignore comments
						}
						else {
							// Oops! Something's wrong
							throw new IllegalStateException("Formatting error. Unrecognized line:\n"+line);
						}
					}

					line = ripc.readLine();
				}
			} catch (IOException e) {
				throw new NoSuchElementException("IOException occured");
			}

			MultipleAlignment m = new MultipleAlignment(RIPCParser.getPDBNames(nextLabels),residues);
			nextLabels = null;
			return m;
		}

		/**
		 * Skips the first level of comments and initializes nextLabels to the first alignment.
		 * @throws IOException 
		 * @throws IOException 
		 */
		private void skipComments() throws IOException {
			String line;
			line = ripc.readLine();

			while( line!=null ) {
				line = line.trim();
				Matcher labels=labelRegex.matcher(line);
				if(labels.matches()) {
					nextLabels = new String[2];
					nextLabels[0] = labels.group(1);
					nextLabels[1] = labels.group(2);
					return;
				}
				Matcher comments=commentRegex.matcher(line);
				if(!comments.matches()) {
					// Oops! Something's wrong
					throw new IllegalStateException("Formatting error. Expected comment or label, found:\n"+line);
				}

				line = ripc.readLine();
			}

			throw new NoSuchElementException();
		}
		

		@Override
		public boolean hasNext()
		{
			return nextLabels != null;
		}



		/**
		 * Not implemented.
		 */
		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}

}
