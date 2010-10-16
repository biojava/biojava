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
package org.biojavax.bio.phylo.io.nexus;

import java.util.List;

import org.biojava.bio.seq.io.ParseException;

/**
 * Builds Nexus characters blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class CharactersBlockBuilder extends NexusBlockBuilder.Abstract
		implements CharactersBlockListener {

	private CharactersBlock block;

	protected void addComment(final NexusComment comment) {
		this.block.addComment(comment);
	}
	
	protected CharactersBlock makeNewBlock() {
		return new CharactersBlock();
	}

	protected NexusBlock startBlockObject() {
		this.block = this.makeNewBlock();
		this.resetStatus();
		return this.block;
	}

	/**
	 * Allowed to be called by DATA subclass.
	 */
	protected void resetStatus() {
		// Nothing to do.
	}

	public void endBlock() {
		// Don't care.
	}

	public void endTokenGroup() {
		// Nothing to do.
	}

	public void addCharLabel(String charLabel) {
		this.block.addCharLabel(charLabel);
	}

	public void addCharState(String charState) {
		this.block.addCharState(charState);
	}

	public void addCharStateKeyword(String charState, String keyword) {
		this.block.addCharStateKeyword(charState, keyword);
	}

	public void addEquate(String symbol, List symbols) {
		this.block.addEquate(symbol, symbols);
	}

	public void addItem(String item) {
		this.block.addItem(item);
	}

	public void addMatrixEntry(String taxa) {
		this.block.addMatrixEntry(taxa);
	}

	public void addState(String state) {
		this.block.addState(state);
	}

	public void addStateLabel(String state, String label) {
		this.block.addStateLabel(state, label);
	}

	public void addSymbol(String symbol) {
		this.block.addSymbol(symbol);
	}

	public void addTaxLabel(String taxLabel) throws ParseException {
		this.block.addTaxLabel(taxLabel);
	}

	public void appendMatrixData(String taxa, Object data) {
		this.block.appendMatrixData(taxa, data);
	}

	public void setCharStateLabel(String charState, String label) {
		this.block.setCharStateLabel(charState, label);
	}

	public void setDataType(String dataType) {
		this.block.setDataType(dataType);
	}

	public void setDimensionsNChar(int dimensionsNChar) {
		this.block.setDimensionsNChar(dimensionsNChar);
	}

	public void setDimensionsNTax(int dimensionsNTax) {
		this.block.setDimensionsNTax(dimensionsNTax);
	}

	public void setEliminateEnd(int eliminateEnd) {
		this.block.setEliminateEnd(eliminateEnd);
	}

	public void setEliminateStart(int eliminateStart) {
		this.block.setEliminateStart(eliminateStart);
	}

	public void setGap(String gap) {
		this.block.setGap(gap);
	}

	public void setInterleaved(boolean interleaved) {
		this.block.setInterleaved(interleaved);
	}

	public void setLabels(boolean labels) {
		this.block.setLabels(labels);
	}

	public void setMatchChar(String matchChar) {
		this.block.setMatchChar(matchChar);
	}

	public void setMissing(String missing) {
		this.block.setMissing(missing);
	}

	public void setRespectCase(boolean respectCase) {
		this.block.setRespectCase(respectCase);
	}

	public void setStatesFormat(String statesFormat) {
		this.block.setStatesFormat(statesFormat);
	}

	public void setTokens(boolean tokens) {
		this.block.setTokens(tokens);
	}

	public void setTransposed(boolean transposed) {
		this.block.setTransposed(transposed);
	}

}
