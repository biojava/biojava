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
package org.biojava.nbio.structure.symmetry.misc;

import org.biojava.nbio.structure.symmetry.utils.BlastClustReader;

import java.util.*;
import java.util.Map.Entry;

public class ProteinComplexSignature {
	private BlastClustReader blastClust = null;
	private String pdbId = "";
	private List<String> chainIds = null;
	private List<ChainSignature> chainSignatures = new ArrayList<ChainSignature>();


	public ProteinComplexSignature(String pdbId, List<String> chainIds, BlastClustReader blastClust) {
		this.pdbId = pdbId;
		this.chainIds = chainIds;
		this.blastClust = blastClust;

		getChainSignatures();
	}

	public String getComplexSignature() {
		StringBuilder builder = new StringBuilder();
		for (ChainSignature s: chainSignatures) {
			builder.append(s.toString());
		}
		return builder.toString();
	}

	public String getCompositionId(String chainId) {
		for (ChainSignature s: chainSignatures) {
			if (s.getChainIds().contains(chainId)) {
				return s.getCompositionId();
			}
		}
		return "";
	}

	public String getComplexStoichiometry() {
		StringBuilder s = new StringBuilder();
		for (ChainSignature c: chainSignatures) {
			s.append(c.getCompositionId());
			if (c.getChainIds().size() >1) {
				s.append(c.getChainIds().size());
			}
		}
		return s.toString();
	}

	public int getSubunitTypeCount() {
		return chainSignatures.size();
	}

	private List<ChainSignature> getChainSignatures() {
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

		Map<String,Integer> mapCounts = new TreeMap<String,Integer>();
		Map<String,List<String>> mapChainIds = new TreeMap<String, List<String>>();

		for (String chainId: chainIds) {
			String rep = blastClust.getRepresentativeChain(pdbId, chainId);
			Integer value = mapCounts.get(rep);
			if (value == null) {
				mapCounts.put(rep, 1);
				List<String> list = new ArrayList<String>();
				list.add(chainId);
				mapChainIds.put(rep, list);
			} else {
				value+=1;
				mapCounts.put(rep, value);
				List<String> list = mapChainIds.get(rep);
				list.add(chainId);
			}
		}


		for (Entry<String, Integer> entry: mapCounts.entrySet()) {
			List<String> chainIds = mapChainIds.get(entry.getKey());
			ChainSignature chainSignature = new ChainSignature(entry.getKey(), entry.getValue(), chainIds);
			chainSignatures.add(chainSignature);
		}

		Collections.sort(chainSignatures);
		for (int i = 0; i < chainSignatures.size(); i++) {
			ChainSignature c = chainSignatures.get(i);
			if (i < alpha.length()) {
				c.setCompositionId(alpha.substring(i,i+1));
			} else {
				c.setCompositionId("?");
			}
		}

		return chainSignatures;
	}

}
