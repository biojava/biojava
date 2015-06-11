package org.biojava.nbio.structure.align.multiple;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;

/**
 * Utility functions for working with multiple alignments. Methods for: sequence alignment calculation.
 * TODO add the new external methods here when implemented.
 * @author Spencer Bliven
 *
 */
public class MultipleAlignmentTools {
	/**
	 * Gets the sequence strings for all blocks in an alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps between blocks are omitted,
	 * while gaps within blocks are represented by '-'.
	 * @param alignment Input alignment
	 * @return a string for each row in the alignment, giving the 1-letter code
	 *  for each aligned residue. All strings have length <tt>alignment.length()</tt>
	 * @see StructureTools#get1LetterCode(String)
	 */
	public static List<String> getSequencesForBlocks(MultipleAlignment alignment) {
		List<String> alnSequences = new ArrayList<String>();
		try {
			List<Atom[]> atoms = alignment.getEnsemble().getAtomArrays();
			for (int row=0; row<alignment.size(); row++){
				String seq = "";
				for (Block b: alignment.getBlocks()){
					for (int res=0; res<b.length(); res++){
						if (b.getAlignRes().get(row).get(res) != null)
							seq += StructureTools.get1LetterCode(atoms.get(row)[res].getGroup().getPDBName());
						else seq += "-";
					}
				}
				alnSequences.add(seq);
			}
		} catch (StructureAlignmentException e) {
			e.printStackTrace();
		}
		return alnSequences;
	}
}
