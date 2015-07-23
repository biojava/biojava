package org.biojava.nbio.structure.align.xml;

import java.io.IOException;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.ScoresCache;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;

/**
 * Helper methods to convert all the hierarchy levels of a MultipleAlignment
 * into an XML format.
 * <p>
 * To convert a MultipleAlignment to an XML String use the 
 * {@link MultipleAlignmentWriter#toXML(MultipleAlignmentEnsemble)} method.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class MultipleAlignmentXMLConverter {
	
	public synchronized static void printXMLensemble(PrettyXMLWriter xml,
			MultipleAlignmentEnsemble ensemble)	throws IOException {

		xml.openTag("MultipleAlignmentEnsemble");
		
		printXMLheader(xml,ensemble);

		for (MultipleAlignment msa:ensemble.getMultipleAlignments()){
			printXMLalignment(xml,msa);
		}
		printXMLscoresCache(xml,ensemble);
		
		xml.closeTag("MultipleAlignmentEnsemble");
	}

	public synchronized static void printXMLalignment(PrettyXMLWriter xml, 
			MultipleAlignment msa) throws IOException {

		xml.openTag("MultipleAlignment");
		
		for(BlockSet bs:msa.getBlockSets()) {
			printXMLblockSet(xml, bs);
		}
		printXMLscoresCache(xml,msa);

		xml.closeTag("MultipleAlignment");
	}

	public synchronized static void printXMLblockSet(PrettyXMLWriter xml,
			BlockSet bs) throws IOException {

		xml.openTag("BlockSet");
		
		for(Block b:bs.getBlocks()) {
			printXMLblock(xml, b);
		}
		
		if (bs.getTransformations() != null){
			for(Matrix4d t:bs.getTransformations()){
				printXMLmatrix4d(xml, t);
			}
		}
		printXMLscoresCache(xml,bs);

		xml.closeTag("BlockSet");
	}

	public synchronized static void printXMLblock(PrettyXMLWriter xml,
			Block b) throws IOException {

		xml.openTag("Block");
		List<List<Integer>> alignment = b.getAlignRes();
		
		for (int pos=0;pos<alignment.get(0).size(); pos++){

			xml.openTag("eqr"+pos);
			for (int str=0; str<alignment.size(); str++){
				xml.attribute("str"+(str+1),alignment.get(str).get(pos)+"");
			}
			xml.closeTag("eqr"+pos);
		}
		printXMLscoresCache(xml,b);
		
		xml.closeTag("Block");
	}

	public synchronized static void printXMLmatrix4d(PrettyXMLWriter xml,
			Matrix4d transform) throws IOException {

		if (transform == null) return;
		xml.openTag("Matrix4d");

		for (int x=0;x<4;x++){
			for (int y=0;y<4;y++){
				String key = "mat"+(x+1)+(y+1);
				String value = transform.getElement(x,y)+"";
				xml.attribute(key,value);
			}
		}
		xml.closeTag("Matrix4d");
	}
	
	public synchronized static void printXMLscoresCache(PrettyXMLWriter xml,
			ScoresCache cache) throws IOException {

		if (cache == null) return;
		xml.openTag("ScoresCache");

		//We need a new tag for every score, we don't know their names
		for (String score:cache.getScores()){
			xml.openTag(score);
			String value = cache.getScore(score)+"";
			xml.attribute("value", value);
			xml.closeTag(score);
		}
		xml.closeTag("ScoresCache");
	}

	public synchronized static void printXMLheader(PrettyXMLWriter xml,
			MultipleAlignmentEnsemble ensemble) throws IOException{

		//Creation properties
		xml.attribute("Algorithm", ensemble.getAlgorithmName());
		xml.attribute("Version", ensemble.getVersion());
		xml.attribute("IOTime", ensemble.getIoTime()+"");
		xml.attribute("CalculationTime", ensemble.getCalculationTime()+"");

		//Structure Identifiers
		xml.openTag("Structures");
		for (int i=0; i<ensemble.size(); i++){
			String name = ensemble.getStructureNames().get(i);
			xml.attribute("name"+(i+1), name);
		}
		xml.closeTag("Structures");
	}
}
