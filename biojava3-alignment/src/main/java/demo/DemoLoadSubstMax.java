package demo;

import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

public class DemoLoadSubstMax {

	public static void main (String[] args){


		// that's the PAM250 matrix (named a bit unclear in AAindex...)
		String max2="DAYM780301";
		SubstitutionMatrix<AminoAcidCompound> substMax2 = SubstitutionMatrixHelper.getMatrixFromAAINDEX(max2);
		System.out.println(substMax2);

		// and here BLOSUM62...
		String max3="HENS920102";
		SubstitutionMatrix<AminoAcidCompound> substMax3 = SubstitutionMatrixHelper.getMatrixFromAAINDEX(max3);
		System.out.println(substMax3);
		
		// This one I developed a while ago to  be optimised for the alignment of distantly related sequences
		String matrixName4 = "PRLA000101";
		SubstitutionMatrix<AminoAcidCompound> substMax4 = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName4);
		System.out.println(substMax4);
		
	}
}
