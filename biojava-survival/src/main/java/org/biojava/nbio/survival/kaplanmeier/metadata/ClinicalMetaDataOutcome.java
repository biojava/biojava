/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.survival.kaplanmeier.metadata;

import org.biojava.nbio.survival.data.WorkSheet;

import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ClinicalMetaDataOutcome {

	/**
	 *
	 * @param worksheet
	 * @param sensorMapColumn
	 * @param censorMap
	 * @param timeColumn
	 * @param timeScale
	 * @param metaDataInfoList
	 * @throws Exception
	 */
	static public void process(WorkSheet worksheet, String sensorMapColumn, LinkedHashMap<String, String> censorMap, String timeColumn, Double timeScale, ArrayList<MetaDataInfo> metaDataInfoList) throws Exception {
		for (MetaDataInfo metaDataInfo : metaDataInfoList) {
			if (metaDataInfo.numeric) {
				metaDataInfo.discreteQuantizer.process(worksheet, metaDataInfo.column);
			}
			metaDataInfo.setDiscreteValues(worksheet);
		}

		for (MetaDataInfo metaDataInfo : metaDataInfoList) {
			int numberValues = metaDataInfo.getNumberDiscreteValues();
			for(int i = 0; i < numberValues; i++){

			}

		}

	}

	/**
	 *
	 * @param args
	 */
	public static void main(String[] args) {

		try {
			LinkedHashMap<String, String> censorMap = new LinkedHashMap<String, String>();
			censorMap.put("a", "0");
			censorMap.put("d", "1");
			censorMap.put("d-d.s.", "1");
			censorMap.put("d-o.c.", "1");
			String timeColumn = "TIME";
			String sensorMapColumn = "last_follow_up_status"; // "survstat3";
			double timeScale = 1.0;
			ArrayList<MetaDataInfo> metaDataInfoList = new ArrayList<MetaDataInfo>();
			metaDataInfoList.add(new MetaDataInfo("age_at_diagnosis", true, new MeanQuantizer()));
			metaDataInfoList.add(new MetaDataInfo("size", true, new MeanQuantizer()));
			metaDataInfoList.add(new MetaDataInfo("lymph_nodes_positive", true, new MeanQuantizer()));
			metaDataInfoList.add(new MetaDataInfo("lymph_nodes_removed", true, new MeanQuantizer()));
			metaDataInfoList.add(new MetaDataInfo("NPI", true, new MeanQuantizer()));
			metaDataInfoList.add(new MetaDataInfo("menopausal_status_inferred"));
			metaDataInfoList.add(new MetaDataInfo("group"));
			metaDataInfoList.add(new MetaDataInfo("grade"));
			metaDataInfoList.add(new MetaDataInfo("stage"));
			metaDataInfoList.add(new MetaDataInfo("ER_IHC_status"));
			metaDataInfoList.add(new MetaDataInfo("HER2_IHC_status"));
			metaDataInfoList.add(new MetaDataInfo("HER2_SNP6_state"));
			metaDataInfoList.add(new MetaDataInfo("cellularity"));
			metaDataInfoList.add(new MetaDataInfo("P53_mutation_status"));
			metaDataInfoList.add(new MetaDataInfo("P53_mutation_type"));
			metaDataInfoList.add(new MetaDataInfo("Pam50Subtype"));
			metaDataInfoList.add(new MetaDataInfo("Genefu"));

			WorkSheet worksheet = WorkSheet.readCSV("/Users/Scooter/scripps/ngs/DataSets/METABRIC/EGAD00010000210/table_S2_revised.txt", '\t');

			ClinicalMetaDataOutcome.process(worksheet, sensorMapColumn, censorMap, timeColumn, timeScale, metaDataInfoList);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
