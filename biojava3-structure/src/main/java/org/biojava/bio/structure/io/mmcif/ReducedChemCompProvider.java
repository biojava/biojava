package org.biojava.bio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** Unlike the {@link DownloadChemCompProvider}, this  {@link ChemCompProvider} does not download any chem comp definitions. 
 * It has access to a limited set of files that are part of the biojava distribution. 
 *   
 * @author Andreas Prlic
 * @since 3.0
 */
public class ReducedChemCompProvider implements ChemCompProvider {
	
	private static final Logger logger = LoggerFactory.getLogger(ReducedChemCompProvider.class);

	public ReducedChemCompProvider(){
		logger.debug("Initialising ReducedChemCompProvider");
	}
	
	
	public ChemComp getEmptyChemComp(){
		ChemComp comp = new ChemComp();
		
		comp.setOne_letter_code("?");
		comp.setPolymerType(PolymerType.unknown);
		comp.setResidueType(ResidueType.atomn);
		return comp;
	}
	
	public ChemComp getChemComp(String recordName) {
		String name = recordName.toUpperCase().trim();
		InputStream inStream = this.getClass().getResourceAsStream("/chemcomp/"+name + ".cif.gz");

		if ( inStream == null){
			//System.out.println("Could not find chem comp: " + name + " ... using generic Chem Comp");
			// could not find the chem comp definition for this in the jar file
			logger.debug("Getting empty chem comp for {}",name);
			ChemComp cc = getEmptyChemComp();
			cc.setId(name);
			return cc;
		}

		MMcifParser parser = new SimpleMMcifParser();

		ChemCompConsumer consumer = new ChemCompConsumer();

		// The Consumer builds up the BioJava - structure object.
		// you could also hook in your own and build up you own data model.
		parser.addMMcifConsumer(consumer);

		try {
			parser.parse(new BufferedReader(new InputStreamReader(new GZIPInputStream(inStream))));

			ChemicalComponentDictionary dict = consumer.getDictionary();

			ChemComp chemComp = dict.getChemComp(name);
			
			return chemComp;
		} catch (IOException e){
			logger.error("IOException caught while reading chem comp {}. Error: {}",name,e.getMessage());
			//e.printStackTrace();

		}
		logger.warn("Problem when loading chem comp {}, will use an empty chem comp for it", name);
		ChemComp cc = getEmptyChemComp();
		cc.setId(name);
		return cc;
	}


}
