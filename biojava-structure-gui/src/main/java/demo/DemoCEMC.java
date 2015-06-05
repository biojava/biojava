package demo;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.cemc.CeMcMain;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Demo for running a CEMC Multiple Structure Alignment and visualizing the results.
 * 
 * @author Aleix Lafita
 * 
 */
public class DemoCEMC {

	public static void main(String[] args) throws IOException, StructureException, StructureAlignmentException, InterruptedException, ExecutionException {
		
		//Specify the structures to align
		//ASP-proteinases (CEMC paper)
		//List<String> names = Arrays.asList("3app", "4ape", "2apr", "5pep", "1psn", "4cms", "1bbs.A", "1smr.A", "2jxr.A", "1mpp", "2asi", "1am5");
		//Protein Kinases (CEMC paper)
		//List<String> names = Arrays.asList("1cdk.A", "1cja.A", "1csn", "1b6c.B", "1ir3.A", "1fgk.A", "1byg.A", "1hck", "1blx.A", "3erk", "1bmk.A", "1kob.A", "1tki.A", "1phk", "1a06");
		//DHFR (Gerstein 1998 paper)
		List<String> names = Arrays.asList("d1dhfa_", "8dfr", "d4dfra_", "3dfr");
		//TIM barrels (MUSTA paper)
		//List<String> names = Arrays.asList("1tim.A", "1vzw", "1nsj", "3tha.A", "4enl", "2mnr", "7tim.A", "1tml", "1btc", "a1piia1", "6xia", "5rub.A", "2taa.B");
		//Helix-bundle (MUSTA paper)
		//List<String> names = Arrays.asList("1bbh.A", "1aep", "1bge.B", "256b.A", "2ccy.A", "2hmz.A", "3ink.C");
		//Calcium Binding (MUSTA paper)
		//List<String> names = Arrays.asList("4cpv", "2scp.A", "2sas", "1top", "1scm.B", "3icb");
		//Serine Rich Proteins SERP (MUSTA paper)
		//List<String> names = Arrays.asList("7api.A", "8api.A", "1hle.A", "1ova.A", "2ach.A", "9api.A", "1psi", "1atu", "1kct", "1ath.A", "1att.A");
		//Serine Proteases (MUSTA paper)
		//List<String> names = Arrays.asList("1cse.E", "1sbn.E", "1pek.E", "3prk", "3tec.E");
		//GPCRs
		//List<String> names = Arrays.asList("4xt3", "4or2.A", "3odu.A", "2z73.A", "4ug2.A");
		//Immunoglobulins (MAMMOTH paper)
		//List<String> names = Arrays.asList("2hla.B", "3hla.B", "1cd8", "2rhe", "1tlk", "1ten", "1ttf");
		//Globins (MAMMOTH and MUSTA papers)
		//List<String> names = Arrays.asList("1mbc", "1hlb", "1thb.A", "1ith.A", "1idr.A", "1dlw", "1kr7.A", "1ew6.A", "1it2.A", "1eco", "3sdh.A", "1cg5.B", "1fhj.B", "1ird.A", "1mba", "2gdm", "1b0b", "1h97.A", "1ash.A", "1jl7.A");
		//Rossman-Fold (POSA paper)
		//List<String> names = Arrays.asList("d1heta2", "d1ek6a_", "d1obfo1", "2cmd", "d1np3a2", "d1bgva1", "d1id1a_", "d1id1a_", "d1oi7a1");
		
		//Load the CA atoms of the structures
		AtomCache cache = new AtomCache();
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (String name:names)	{
			atomArrays.add(cache.getAtoms(name));
		}
		
		//Here the multiple structural alignment algorithm comes in place to generate the alignment object
		CeMcMain algorithm = new CeMcMain();
		MultipleAlignment result = algorithm.align(atomArrays);
		result.getEnsemble().setStructureNames(names);
		
		//Information about the alignment
		result.getEnsemble().setAlgorithmName(algorithm.getAlgorithmName());
		result.getEnsemble().setVersion(algorithm.getVersion());
        
		StructureAlignmentDisplay.display(result);
	}
}
