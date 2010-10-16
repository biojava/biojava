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
 * Parses files formatted like the CPDB dataset.
 * <i>Lo et al. CPDB: a database of circular permutation in proteins. Nucleic Acids Res (2009) vol. 37 (Database issue) pp. D328-32</i>
 * <br/>Available through the <a href="http://sarst.life.nthu.edu.tw/cpdb/">CPDB</a>
 * or at src/test/resources/align/benchmarks/CPDB*.txt.
 * Machine readable files courtesy of Wei-Cheng Lo and Ping-Chiang Lyu
 * <p>
 * See {@link CDCPIterator} for info about the format.
 * @author Spencer Bliven
 *
 */
public class CPDBParser implements MultipleAlignmentParser
{
	private String filename;
	public CPDBParser(String filename) {
		this.filename=filename;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		String filename = "src/test/resources/align/benchmarks/CPDB_CPpairsAlignments_withCPSiteRefinement.txt";
		CPDBParser parser = new CPDBParser(filename);
		try {
		for(MultipleAlignment ma : parser) {
			System.out.println(ma.display());
		}
		} catch(IllegalStateException e) {
			e.printStackTrace();
			System.exit(1);
			return;
		}
	}

	//@Override
	public Iterator<MultipleAlignment> iterator()
	{
		try {
			return new CPDBIterator(filename);
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}


	/**
	 * Parses a CPDB alignment file.
	 * <p>
	 * Format:<br/>
	 * Comments start with ';'.
	 * Each alignment begins with a line starting with '#', then contains a
	 * header section of the form 'name : value'. The 'Qry' and 'Sbj' fields
	 * give the PDBID and chain for the alignment. Finally, the alignment pairs
	 * are listed in the 'AliResPairs' section.
	 * <p>
	 * Example:<br/><pre>
	 * # 1	121p__cp137	1htwA
	 * Qry        : 121p_ (166 a.a.)
	 * Sbj        : 1htwA (158 a.a.)
	 * CPSite@Qry : 137
	 * RMSD       : 4.6526
	 * Str. Div.  : 12.6879
	 * CP Score   : 0.2427
	 * Ali. Size  : 83 a.a.
	 * Identity   : 8/158 (5.1%)
	 * 
	 * Qry AliRes : YGIPYIETSAKT-RQGVEDAFYTLVREIRQH--MTEYKLVVVGAGGVGKSALTIQLI...
	 * Sbj AliRes : MESLTQYIP---DEFSMLRFGKKFAEILLKLHTEKAIMVYLNGDLGAGKTTLTRGML...
	 *
	 * AliResPairs:
	 * ; Qry Res <-> Sbj Res
	 *   137 TYR Y M MET 1             
	 *   138 GLY G E GLU 2  
	 * </pre>           
	 *
	 * @author Spencer Bliven
	 *
	 */
	protected static class CPDBIterator implements Iterator<MultipleAlignment> {
		private BufferedReader cpdb;
		private String[] nextLabels;
		private String[] nextChains;

		private static final Pattern headerRegex =
			Pattern.compile("^(.*?)\\s*:\\s*(.*?)$");
		private static final Pattern pdbIdentifierRegex =
			Pattern.compile("(....)(.)\\s+\\(.*\\)"); //for Qry/Sbj headers
		private static final Pattern pairRegex = 
			Pattern.compile("^(-?\\d+\\S?)\\s+([A-Z]{3}|\\?)\\s+[A-Z?]\\s+" +
					"[A-Z?]\\s+([A-Z]{3}|\\?)\\s+(-?\\d+\\S?)$");
		private static final Pattern commentRegex = 
			Pattern.compile("^(?:;.*|$)");
		private static final Pattern identifierRegex = 
			Pattern.compile("^#.*$");

		public CPDBIterator(String filename) throws IOException {
			cpdb = new BufferedReader(new FileReader(filename));
			assert(cpdb.markSupported());
			nextLabels = null;
			nextChains = null;
			skipComments();
		}

		private static String formatPDBID(String pdbID, String chain) {
			if(chain.equals("_")) {
				return pdbID;
			} else {
				return pdbID+"."+chain;
			}
		}
		
		//@Override
		public MultipleAlignment next()
		{
			List<List<PDBResidue>> residues = new ArrayList<List<PDBResidue>>(2);
			residues.add(new LinkedList<PDBResidue>());
			residues.add(new LinkedList<PDBResidue>());


			String line;
			try {
				line = cpdb.readLine();

				while( line!=null ) {
					line = line.trim();
					Matcher header = headerRegex.matcher(line);
					if(header.matches()) {
						// Store Qry/Sbj PDB IDs
						String name = header.group(1);
						String label = header.group(2);

						if(name.equals("Qry")) {
							Matcher pdbIDs = pdbIdentifierRegex.matcher(label);
							if(! pdbIDs.matches()) {
								throw new IllegalStateException("Formatting error. Expected a pdbID in '"+label+"'");
							}
							
							this.nextLabels[0] = formatPDBID(pdbIDs.group(1),pdbIDs.group(2));
							this.nextChains[0] = pdbIDs.group(2);
						}
						else if(name.equals("Sbj")) {
							Matcher pdbIDs = pdbIdentifierRegex.matcher(label);
							if(! pdbIDs.matches()) {
								throw new IllegalStateException("Formatting error. Expected a pdbID in '"+label+"'");
							}
							this.nextLabels[1] = formatPDBID(pdbIDs.group(1),pdbIDs.group(2));
							this.nextChains[1] = pdbIDs.group(2);
						}
						
					} else {
						Matcher pair = pairRegex.matcher(line);
						if(pair.matches()) {
							String aa1 = pair.group(2);
							String aa2 = pair.group(3);
							
							if(aa1.equals("?")) {
								aa1 = null;
							}
							if(aa2.equals("?")) {
								aa2 = null;
							}

							String pdb1 = pair.group(1);
							String pdb2 = pair.group(4);
							
							residues.get(0).add(new PDBResidue(pdb1,nextChains[0],aa1));
							residues.get(1).add(new PDBResidue(pdb2,nextChains[1],aa2));
						}
						else {
							Matcher identifier=identifierRegex.matcher(line);
							if(identifier.matches()) {
								if(nextLabels[0]==null) {
									throw new IllegalStateException("Found a new record before a Qry line: "+line);
								}
								if(nextLabels[1]==null) {
									throw new IllegalStateException("Found a new record before a Sbj line: "+line);
								}
								MultipleAlignment m = new MultipleAlignment(nextLabels,residues);

								nextLabels = new String[] { null, null };
								nextChains = new String[] { null, null };

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
					}

					line = cpdb.readLine();
				}
				
			} catch (IOException e) {
				throw new NoSuchElementException("IOException occured");
			}

			MultipleAlignment m = new MultipleAlignment(nextLabels,residues);
			nextLabels = null;
			nextChains = null;
			return m;
		}

		/**
		 * Skips the first level of comments and initializes nextLabels.
		 * @throws IOException 
		 * @throws IOException 
		 */
		private void skipComments() throws IOException {
			String line;
			line = cpdb.readLine();

			while( line!=null ) {
				line = line.trim();
				Matcher identifier=identifierRegex.matcher(line);
				if(identifier.matches()) {
					nextLabels = new String[] { null, null };
					nextChains = new String[] { null, null };
					return;
				}
				Matcher comments=commentRegex.matcher(line);
				if(!comments.matches()) {
					// Oops! Something's wrong
					throw new IllegalStateException("Formatting error. Expected comment or header, found:\n"+line);
				}

				line = cpdb.readLine();
			}

			throw new NoSuchElementException();
		}


		//@Override
		public boolean hasNext()
		{
			return nextLabels != null;
		}



		/**
		 * Not implemented.
		 */
		//@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}

}
