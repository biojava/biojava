package org.biojava.nbio.structure.align.multiple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Utility functions for working with {@link MultipleAlignment}. 
 * <p>
 * Supported functions:
 * <ul><li>Multiple Sequence Alignment Calculation Methods
 * <li>Map from Sequence Alignment Position to Structure Atom
 * <li>Map from Sequence Alignment Position to Block Number
 * </ul>
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentTools {
	
	/**
	 * Calculate the sequence alignment Strings for the whole alignment. This method creates a sequence 
	 * alignment where aligned residues are in uppercase and unaligned residues are in lowercase, thus
	 * providing a more compact way to represent the alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps are represented by '-'. Separation between different 
	 * Blocks is indicated by a gap in all positions, meaning that there is a possible discontinuity.
	 * 
	 * @param alignment input MultipleAlignment
	 * @param mapSeqToStruct provides a link from the sequence alignment position to the structure alignment 
	 * 		  position. Specially designed for the GUI. Has to be initialized previously and will be overwritten.
	 * @return a string for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 */
	public static List<String> getSequenceAlignment(MultipleAlignment alignment, List<Integer> mapSeqToStruct) {

		//Initialize sequence variables
		List<String> alnSequences = new ArrayList<String>();
		for (int str=0; str<alignment.size(); str++) alnSequences.add("");
		mapSeqToStruct.clear();
		List<Atom[]> atoms = alignment.getEnsemble().getAtomArrays();
		int globalPos = -1;
				
		//Initialize helper variables in constucting the sequence alignment: freePool and blockStarts
		List<SortedSet<Integer>> freePool = new ArrayList<SortedSet<Integer>>();
		List<SortedSet<Integer>> blockStarts = new ArrayList<SortedSet<Integer>>();
		List<List<Integer>> aligned = new ArrayList<List<Integer>>();
		
		//Generate freePool residues from the ones not aligned
		for (int i=0; i<alignment.size(); i++){
			List<Integer> residues = new ArrayList<Integer>();
			freePool.add(new TreeSet<Integer>());
			blockStarts.add(new TreeSet<Integer>());
			for (BlockSet bs : alignment.getBlockSets()){
				for (Block b : bs.getBlocks()){
					boolean first = true;
					for (int l=0; l<b.length(); l++){
						Integer residue = b.getAlignRes().get(i).get(l);
						if (residue != null){
							if (first) blockStarts.get(i).add(residue);
							residues.add(residue);
							first = false;
						}
					}
				}
			}
			aligned.add(residues);
		}
		//Add any residue not aligned to the free pool for every structure
		for (int i=0; i<alignment.size(); i++){
			for (int k=0; k<atoms.get(i).length; k++){
				if (!aligned.get(i).contains(k)) freePool.get(i).add(k);
			}
		}
		
		//Loop through all the alignment Blocks in the order given
		for (int b=0; b<alignment.getBlocks().size(); b++){
			if (b!=0){
				//Add a gap to all structures in order to separate visually the blocks in the alignment
				for (int str=0; str<alignment.size(); str++) alnSequences.set(str,alnSequences.get(str).concat("-"));
				mapSeqToStruct.add(-1); //means no aligned position
			}
			//Store the previous position added to the sequence alignment for this structure
			int[] previousPos = new int[alignment.size()];
			Arrays.fill(previousPos, -1);
			//Store provisional characters
			char[] provisionalChar = new char[alignment.size()];
			Arrays.fill(provisionalChar, '-');
			//Loop through all the alignment positions in the Block
			for (int pos=0; pos<alignment.getBlocks().get(b).length(); pos++){
				globalPos++;
				boolean gaps = true;  //If any structure is not consecutive with the previousPos
				//While the next position cannot be considered because there are still non consecutive residues
				while (gaps){
					gaps = false;
					//Loop through all the structures
					for (int str=0; str<alignment.size(); str++){
						//If it is the first position or before it was null
						if (previousPos[str] == -1){
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) provisionalChar[str] = '-';
							else provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());
						}
						else {
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) {
								if (freePool.get(str).contains(previousPos[str]+1))
									provisionalChar[str] = Character.toLowerCase(StructureTools.get1LetterCode(atoms.get(str)[previousPos[str]+1].getGroup().getPDBName()));
								else provisionalChar[str] = '-';
							}
							else if (previousPos[str]+1 == residue){
								provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());
							} else{
								provisionalChar[str] = ' ';  //This means there is a spacing (non-consecutive)
								gaps = true;
							}
						}
					}//End all structure analysis
					if (gaps){
						for (int str=0; str<alignment.size(); str++){
							if (provisionalChar[str] == ' ') {  //It means this residue was the non-consecutive one
								alnSequences.set(str,alnSequences.get(str).concat(""+Character.toLowerCase(StructureTools.get1LetterCode(atoms.get(str)[previousPos[str]+1].getGroup().getPDBName()))) );
								previousPos[str] ++;
							} else { //insert a gap otherwise and do not change the index
								alnSequences.set(str,alnSequences.get(str).concat("-"));
							}
						}
						mapSeqToStruct.add(-1); //meaning that this is an unaligned position
					} 
					else {  //Append the provisional and update the indices otherwise
						for (int str=0; str<alignment.size(); str++){
							alnSequences.set(str,alnSequences.get(str).concat(""+provisionalChar[str]));
							if (provisionalChar[str] != '-') {
								if (alignment.getBlocks().get(b).getAlignRes().get(str).get(pos)==null) previousPos[str]++;
								else previousPos[str] = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							}
						}
						mapSeqToStruct.add(globalPos);
					}
				}
			}
			//Now add the unaligned residues in between Blocks (lowercase), using intervals
			int[] blockEnds = new int[alignment.size()];
			for (int str=0; str<alignment.size(); str++){
				for (int res:blockStarts.get(str)){
					if (previousPos[str] > res) blockEnds[str] = res;
					else {
						blockEnds[str] = res;
						break;
					}
				}
			}
			//Now we have the ending residues of each Block before the others, put the residues sequentially
			boolean allGaps = false;
			while (!allGaps){
				allGaps = true;
				for (int str=0; str<alignment.size(); str++){
					if (previousPos[str]+1 < blockEnds[str]){
						alnSequences.set(str,alnSequences.get(str).concat(""+Character.toLowerCase(StructureTools.get1LetterCode(atoms.get(str)[previousPos[str]+1].getGroup().getPDBName()))) );
						previousPos[str]++;
						allGaps = false;
					} else alnSequences.set(str,alnSequences.get(str).concat("-"));
				}
				mapSeqToStruct.add(-1);
			}
		}
		return alnSequences;
	}
	
	/**
	 * Calculate the sequence alignment Strings for the whole alignment. This method creates a sequence 
	 * alignment where aligned residues are in uppercase and unaligned residues are in lowercase, thus
	 * providing a more compact way to represent the alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps are represented by '-'. Separation between different 
	 * Blocks is indicated by a gap in all positions, meaning that there is a possible discontinuity.
	 * 
	 * @param alignment input MultipleAlignment
	 * @return String for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 */
	public static List<String> getSequenceAlignment(MultipleAlignment alignment) {
		return getSequenceAlignment(alignment, new ArrayList<Integer>());
	}
	
	/**
	 * Calculate the sequence alignment Strings for the alignment Blocks in an alignment. This method 
	 * creates a sequence alignment where all residues are in uppercase and a residue alone with gaps
	 * in all the other structures represents unaligned residues. Because of this latter constraint only
	 * the residues within the Blocks are represented, for a more compact alignment. For a sequence alignment of the full protein use
	 * {@link #getSequenceAlignment(MultipleAlignment)}.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps between blocks are omitted,
	 * while gaps within blocks are represented by '-'. Separation between different Blocks is
	 * indicated by a gap in all positions, meaning that there is something unaligned inbetween.
	 *
	 * @param alignment input MultipleAlignment
	 * @param mapSeqToStruct provides a link from the sequence alignment position to the structure alignment 
	 * 		  position. Specially designed for the GUI. Has to be initialized previously and will be overwritten.
	 * @return a string for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 */
	public static List<String> getBlockSequenceAlignment(MultipleAlignment alignment, List<Integer> mapSeqToStruct) {

		//Initialize sequence variables
		List<String> alnSequences = new ArrayList<String>();
		for (int str=0; str<alignment.size(); str++) alnSequences.add("");
		mapSeqToStruct.clear();
		List<Atom[]> atoms = alignment.getEnsemble().getAtomArrays();
		int globalPos = -1;
		
		//Loop through all the alignment Blocks in the order given
		for (int b=0; b<alignment.getBlocks().size(); b++){
			if (b!=0){
				//Add a gap to all structures in order to separate visually the blocks in the alignment
				for (int str=0; str<alignment.size(); str++) alnSequences.set(str,alnSequences.get(str).concat("-"));
				mapSeqToStruct.add(-1); //means no aligned position
			}
			//Store the previous position added to the sequence alignment for this structure
			int[] previousPos = new int[alignment.size()];
			Arrays.fill(previousPos, -1);
			//Store provisional characters
			char[] provisionalChar = new char[alignment.size()];
			Arrays.fill(provisionalChar, '-');
			//Loop through all the alignment positions in the Block
			for (int pos=0; pos<alignment.getBlocks().get(b).length(); pos++){
				globalPos++;
				boolean gaps = true;  //If any structure is not consecutive with the previousPos
				//While the next position cannot be considered because there are still non consecutive residues
				while (gaps){
					gaps = false;
					//Loop through all the structures
					for (int str=0; str<alignment.size(); str++){
						//If it is the first position or before it was null
						if (previousPos[str] == -1){
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) provisionalChar[str] = '-';
							else provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());		
						}
						else{
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) provisionalChar[str] = '-';
							else if (previousPos[str]+1 == residue){
								provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());
							}
							else{
								provisionalChar[str] = ' ';  //This means there is a spacing (non-consecutive)
								gaps = true;
							}
						}
					}//End all structures analysis
					if (gaps){
						for (int str=0; str<alignment.size(); str++){
							if (provisionalChar[str] == ' ') {  //It means this residue was a non-consecutive one, put gaps elsewere
								for (int str2=0; str2<alignment.size(); str2++){
									if (str==str2) 
										alnSequences.set(str2,alnSequences.get(str2).concat(""+StructureTools.get1LetterCode(atoms.get(str2)[previousPos[str]+1].getGroup().getPDBName())));
									else alnSequences.set(str2,alnSequences.get(str2).concat("-")); //insert a gap otherwise and do not change the index
								}
								mapSeqToStruct.add(-1); //meaning that this is an unaligned position
								previousPos[str] += 1;
							}
						}
					}
					else {  //Append the provisional and update the indices otherwise
						for (int str=0; str<alignment.size(); str++){
							alnSequences.set(str,alnSequences.get(str).concat(""+provisionalChar[str]));
							if (provisionalChar[str] != '-') previousPos[str] = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
						}
						mapSeqToStruct.add(globalPos);
					}
				}
			}
		}
		return alnSequences;
	}
	
	/**
	 * Calculate the sequence alignment Strings for the alignment Blocks in an alignment. This method 
	 * creates a sequence alignment where all residues are in uppercase and a residue alone with gaps
	 * in all the other structures represents unaligned residues. Because of this latter constraint only
	 * the residues within the Blocks are represented, for a more compact alignment. For a sequence alignment of the full protein use
	 * {@link #getSequenceAlignment(MultipleAlignment)}.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps between blocks are omitted,
	 * while gaps within blocks are represented by '-'. Separation between different Blocks is
	 * indicated by a gap in all positions, meaning that there is something unaligned inbetween.
	 * 
	 * @param alignment input MultipleAlignment
	 * @return String for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 */
	public static List<String> getBlockSequenceAlignment(MultipleAlignment alignment) {
		return getBlockSequenceAlignment(alignment, new ArrayList<Integer>());
	}
	
   /**
    * Returns the Atom of the specified structure that is aligned in the sequence alignment position specified.
    * 
    * @param multAln the MultipleAlignment object from where the sequence alignment has been generated
    * @param mapSeqToStruct the mapping between sequence and structure generated with the sequence alignment
    * @param structure the structure index of the alignment (row)
    * @param sequencePos the sequence alignment position (column)
    * @return Atom the atom in that position or null if there is a gap
    */
   public static Atom getAtomForSequencePosition(MultipleAlignment multAln, List<Integer> mapSeqToStruct, int structure, int sequencePos) {
	   
	   int seqPos = mapSeqToStruct.get(sequencePos);
	   //Check if the position selected is an aligned position
	   if (seqPos == -1) return null;
	   else {
		   Atom a = null;
		   //Calculate the corresponding structure position (by iterating all Blocks)
		   int sum = 0;
		   for (Block b:multAln.getBlocks()){
			   if (sum+b.length()<=seqPos) {
				   sum += b.length();
				   continue;
			   } else {
				   for (Integer p:b.getAlignRes().get(structure)){
					   if (sum == seqPos) {
						   if (p!= null) a = multAln.getEnsemble().getAtomArrays().get(structure)[p];
						   break;
					   }
					   sum++;
				   }
				   break;
			   }
		   }
		   return a;
	   }
   }
   
   /**
    * Returns the block number of a specified position in the sequence alignment.
    * 
    * @param multAln the MultipleAlignment object from where the sequence alignment has been generated
    * @param mapSeqToStruct the mapping between sequence and structure generated with the sequence alignment
    * @param sequencePos the position in the sequence alignment (column)
    * @return int the block index, or -1 if the position is not aligned
    */
   public static int getBlockForSequencePosition(MultipleAlignment multAln, List<Integer> mapSeqToStruct, int sequencePos){
	   
	   int seqPos = mapSeqToStruct.get(sequencePos);
	   //Check if the position selected is an aligned position
	   if (seqPos == -1) return -1;
	   else {
		   //Calculate the corresponding block (by iterating all Blocks)
		   int sum = 0;
		   int block = 0;
		   for (Block b:multAln.getBlocks()){
			   if (sum+b.length()<=seqPos) {
				   sum += b.length();
				   block++;
				   continue;
			   } else break;
		   }
		   return block;
	   }
   }
   
   /**
    * The average residue distance Matrix contains the average distance from each residue to all 
    * other residues aligned with it. <p>
    * Complexity: T(n,l) = O(l*n^2), if n=number of structures and l=alignment length.
    * 
    * @param alignment MultipleAlignment
    * @return Matrix containing all average residue distances in alignmed columns
    */
   public static Matrix getAverageResidueDistances(MultipleAlignment alignment){
	   
	   
	   
	   
	   
	   
	   return null;
   }
}