package org.biojava.nbio.structure.align.multiple;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.multiple.util.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by andreas on 1/28/16.
 */
public class MultipleAlignmentUtils {

    private static final Logger logger =
            LoggerFactory.getLogger(MultipleAlignmentUtils.class);



    /** New structures are downloaded if they were
     *  not cached in the alignment and they are entirely
     *  transformed here with the superposition information
     *  in the Multiple Alignment.
     *
     * @param multAln
     * @return list of transformed AtomArrays
     * @throws StructureException
     */


    public static List<Atom[]> getRotatedAtoms(MultipleAlignment multAln) throws StructureException{

        int size = multAln.size();

        List<Atom[]> atomArrays = multAln.getAtomArrays();
        for (int i=0; i<size; i++){
            if (atomArrays.get(i).length < 1)
                throw new StructureException(
                        "Length of atoms arrays is too short! Size: "
                                + atomArrays.get(i).length);
        }

        List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();

        //TODO implement independent BlockSet superposition of the structure
        List<Matrix4d> transf = multAln.getBlockSet(0).getTransformations();

        if(transf == null) {

            logger.error("Alignment Transformations are not calculated. "
                    + "Superimposing to first structure as reference.");

            multAln = multAln.clone();
            MultipleSuperimposer imposer = new ReferenceSuperimposer();
            imposer.superimpose(multAln);
            transf = multAln.getBlockSet(0).getTransformations();
            assert(transf != null);
        }

        //Rotate the atom coordinates of all the structures
        for (int i=0; i<size; i++){
            //TODO handle BlockSet-level transformations
            //make sure this method has the same behavior as the other display.
            //-SB 2015-06

            //Assume all atoms are from the same structure
            Structure displayS = atomArrays.get(i)[0].getGroup().
                    getChain().getParent().clone();
            //Get all the atoms and include ligands and hetatoms
            Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
            List<Group> hetatms = StructureTools.getUnalignedGroups(rotCA);
            for (Group g:hetatms){
                rotCA = Arrays.copyOf(rotCA, rotCA.length + 1);
                rotCA[rotCA.length - 1] = g.getAtom(0);
            }

            //Transform the structure to ensure a full rotation in the display
            Calc.transform(displayS, transf.get(i));
            rotatedAtoms.add(rotCA);
        }

    return rotatedAtoms;
    }
}
