package demo;

import org.biojava.bio.structure.io.mmcif.DownloadChemCompProvider;

public class ChemCompDistribution {

	public static void main(String[] args){

		DownloadChemCompProvider c = new DownloadChemCompProvider();
		c.setDownloadAll(true);
		c.checkDoFirstInstall();

	}
}
