/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import java.util.Arrays;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;

/**
 * Stores the alignment between multiple proteins. For a set of proteins p1,p2,...,
 * an alignment is a set of tuples (aa1,aa2,...) specifying a residue from each protein.
 * These residues are defined to be equivalent.
 * Gaps are not allowed; there is no way to represent an alignment between some 
 * proteins which excludes other proteins.
 * @author Spencer Bliven
 *
 */
public class MultipleAlignment
{
	/**
	 * The pdb code: A residue number plus the insertion code
	 * @see Group#getPDBCode()
	 */
	private PDBResidue[][] residues;
	private String[] names;

	/**
	 * Represents a comparison between <i>n</i> proteins.
	 * @param residues An <i>n</i>x<i>m</i> Matrix. 
	 *  Each row should represent a single structure. 
	 *  Each column gives a list of {@link PDBResidue}s from each protein which should be
	 *  aligned.
	 * @param proteinNames An identifier for each protein (eg '1TIM.A')
	 */
	public MultipleAlignment(String[] proteinNames, PDBResidue[][] residues)
	{	
		if(residues.length>0) {
			for(PDBResidue[] prot : residues) {
				if( prot.length != residues[0].length) {
					throw new IllegalArgumentException("Argument must be a rectangular array.");
				}
			}
		}
		this.residues = residues;
		this.names = proteinNames;
	}

	/**
	 * Alternate constructor for 2D Lists instead of arrays
	 * @param proteinNames
	 * @param residues
	 * @see #MultipleAlignment(String[], PDBResidue[][])
	 */
	public MultipleAlignment(String[] proteinNames, List<List<PDBResidue>> residues) {
		this(proteinNames,listToArray(residues));
	}
	/**
	 * Helper function to convert 2D lists into arrays
	 */
	private static PDBResidue[][] listToArray(List<List<PDBResidue>> list) {
		PDBResidue[][] arr =  new PDBResidue[list.size()][];
		int i=0;
		for(List<PDBResidue> res : list) {
			arr[i] = res.toArray(new PDBResidue[res.size()]);
			i++;
		}
		
		return arr;
	}

	/**
	 * @return An <i>n</i>x<i>m</i> Matrix. 
	 *  Each row should represent a single structure. 
	 *  Each column gives a list of {@link PDBResidue}s from each protein which should be
	 *  aligned.
	 */
	public PDBResidue[][] getAlignmentResidues()
	{
		return residues;
	}

	/**
	 * @param residues An <i>n</i>x<i>m</i> Matrix. 
	 *  Each row should represent a single structure. 
	 *  Each column gives a list of {@link PDBResidue}s from each protein which should be
	 *  aligned.
	 */
	public void setAlignmentResidues(PDBResidue[][] residues)
	{
		if(residues.length>0) {
			for(PDBResidue[] prot : residues) {
				if( prot.length != residues[0].length) {
					throw new IllegalArgumentException("Argument must be a rectangular array.");
				}
			}
		}
		this.residues = residues;
	}

	/**
	 * @return the proteinLabels
	 */
	public String[] getNames()
	{
		return names;
	}

	/**
	 * @param names the proteinLabels to set
	 */
	public void setProteinLabels(String[] proteinIDs)
	{
		this.names = proteinIDs;
	}

	/**
	 * Converts the PDB coordinates stored locally (accessible via {@link #getAlignmentMatrix()})
	 * into internal residue numbers.
	 * <p>
	 * The conversion can be run for just a subset of the MultipleAlignment by 
	 * specifying the appropriate protein IDs.
	 * @param structures A list of structures for proteins in this MultipleAlignment.
	 *  Used to map PDB residue numbers to internal coordinates.
	 * @param pIDs For each entry in structures, specifies what alignment
	 *  protein each structure corresponds to. This should be the same length as
	 *  structures, and should be a subset of the values in {@link #getNames()}
	 * @return A matrix with a row for each item in structures. Columns correspond
	 *  to aligned residues between the specified proteins. 
	 * @throws StructureException If a pdb number does not appear in the corresponding Structure.
	 */
	public int[][] getAlignmentMatrix(String[] pIDs, List<Atom[]> structures) throws StructureException {
		int[] indices = new int[pIDs.length];
		for(int i=0;i<pIDs.length;i++) {
			int index;
			for(index=0; index<this.names.length; index++) {
				if( this.names[index].equals( pIDs[i] )) {
					break;
				}
			}
			if(index==this.names.length) {
				throw new IllegalArgumentException("First argument expected to " +
						"be a subset of getProteinIDs(). \""+pIDs[i]+"\" not found.");
			}
			indices[i]=index;
		}
		return getAlignmentMatrix(indices,structures);
	}
	
	/**
	 * Converts the PDB coordinates of all proteins into internal residue numbers
	 * @param structures One structure for each protein
	 * @return
	 * @throws StructureException
	 */
	public int[][] getAlignmentMatrix(List<Atom[]> structures) throws StructureException {
		if(structures.size() != this.residues.length) {
			throw new IllegalArgumentException(
					String.format("Error: %d structures provided, but alignment is length %d\n",
							structures.size(),residues.length) );
		}
		int[] indices = new int[this.names.length];
		for(int i=0;i< indices.length;i++) {
			indices[i]=i;
		}
		return getAlignmentMatrix(indices, structures);
	}

	/**
	 * Converts the PDB coordinates stored locally (accessible via {@link #getAlignmentMatrix()})
	 * into internal residue numbers.
	 * <p>
	 * The conversion can be run for just a subset of the MultipleAlignment by 
	 * specifying the appropriate protein indices.
	 * @param structures A list of structures for proteins in this MultipleAlignment.
	 *  Used to map PDB residue numbers to internal coordinates.
	 * @param proteinIndices For each entry in structures, specifies what alignment
	 *  protein each structure corresponds to. This should be the same length as
	 *  structures, and should specify the row index of each structure into the 
	 *  {@link #getAlignmentResidues()} matrix.
	 * @return A matrix with a row for each item in structures. Columns correspond
	 *  to aligned residues between the specified proteins.
	 * @throws StructureException If a pdb number does not appear in the corresponding Structure.
	 */
	public int[][] getAlignmentMatrix(int[] proteinIndices, List<Atom[]> structures) throws StructureException {
		if(proteinIndices.length != structures.size() ) {
			throw new IllegalArgumentException("Arguments are expected to be the same length");
		}

		int[][] alignMat = new int[proteinIndices.length][];

		for(int prot=0;prot<proteinIndices.length;prot++) {
			int protIndex = proteinIndices[prot];
			alignMat[prot] = new int[residues[protIndex].length];
			Atom[] struct = structures.get(protIndex);

			int atomIndexGuess = 0;//reduce the time to identify Atoms by assuming sequentiality
			for(int res=0;res<residues[protIndex].length;res++) {
				//				Chain chain = struct.getChainByPDB(chains[prot][res]);
				//				Group group = chain.getGroupByPDB(Integer.toString(residues[prot][res]));
				//				alignMat[prot][res] = group.get
				int groupNr = findGroup(struct,residues[protIndex][res],atomIndexGuess);
				// Check the Group number exists
				if(groupNr<0) {
					throw new StructureException(String.format(
							"Unable to locate residue %s in %s.",
							residues[protIndex][res].toString(), struct[0].getParent().getParent().getParent().getPDBCode() ));
				}
				alignMat[prot][res] = groupNr;
			}
		}


		return alignMat;
	}

	/**
	 * Locates the 1st instance of a group with the specified parameters within
	 * an array of atoms and returns the index thereof.
	 * <p>If the amino acid for pdbCode is specified, requires that it match.
	 * <p>Requires O(atoms.length) time.
	 * @param atoms The array of CA atoms to search in. 
	 * @param residueNr The PDB residue number to search for
	 * @param chain The chain to search for
	 * @param atomIndexGuess Position in the array to start the search. Search 
	 * will wrap to the start when the end of the array is reached.
	 * @return
	 */
	private static int findGroup(Atom[] atoms, PDBResidue pdbCode, int atomIndexGuess) {
		int pos=atomIndexGuess;
		do {
			Group g = atoms[pos].getParent();
			Chain c = g.getParent();
			if( ("_".equals(pdbCode.getChain()) || c.getName().equals(pdbCode.getChain()) ) &&
					g.getPDBCode().equals(pdbCode.getResidueCode()) &&
					(pdbCode.getAaName() == null || pdbCode.getAaName().equals(g.getPDBName()) )) {
				return pos;
			}
			pos++;
			pos%=atoms.length;
		} while(pos!=atomIndexGuess);

		return -1;
	}


	/**
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "MultipleAlignment between " + names.length + " proteins: "
		+ Arrays.toString(names);
	}

	public String display() {
		if(residues == null || residues.length<1)
			return "";

		StringBuilder str = new StringBuilder();
		join(str, names, "\t");
		str.append("\n");

		for(int res=0;res<residues[0].length;res++) {
			str.append( residues[0][res].toString());

			for(int prot=1;prot<residues.length;prot++) {
				str.append("\t").append( residues[prot][res].toString());
			}
			str.append("\n");
		}

		return str.toString();
	}

	/**
	 * Joins a list of objects together into a string. The string representation
	 * of each object is calculated from the {@link String#valueOf} method.
	 * 
	 * @param str A stringBuilder object to add the created string to
	 * @param objs A list of objects to join.
	 * @param separator The delimiter string
	 */
	private static void join( StringBuilder str, Object[] objs, String separator )
	{
		if ( objs == null || ( objs.length<1 ) )
			return;
		str.append( objs[0] );
		for(int i=1;i<objs.length;i++) {
			str.append( separator ).append( objs[i] );
		}
	}
	
	public int size() {
		if( residues.length<1 ) {
			return 0;
		}
		return residues[0].length;
	}

}
