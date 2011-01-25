package demo;

import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;

public class ChemCompDistribution {

	public static void main(String[] args){

		
		DownloadChemCompProvider.setPath("/Users/andreas/WORK/PDB/");
		DownloadChemCompProvider c = new DownloadChemCompProvider();
		c.setDownloadAll(true);
		c.checkDoFirstInstall();

	}
}
