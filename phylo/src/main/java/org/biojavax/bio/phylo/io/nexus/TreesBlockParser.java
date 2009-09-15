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

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.TreesBlock.NewickTreeString;

/**
 * Parses Nexus taxa blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class TreesBlockParser extends NexusBlockParser.Abstract {

	public void beginComment() {
		if (!this.expectingTreeRooted)
			super.beginComment();
	}

	public void commentText(String comment) throws ParseException {
		if (this.expectingTreeRooted)
			this.currentTreeRoot.append(comment);
		else
			super.commentText(comment);
	}

	public void endComment() {
		if (!this.expectingTreeRooted)
			super.endComment();
	}

	private boolean expectingTranslate;

	private boolean expectingTaxLabel;

	private boolean expectingTaxName;

	private boolean expectingTree;

	private boolean expectingTreeName;

	private boolean expectingTreeEquals;

	private boolean expectingTreeContent;

	private boolean expectingTreeRooted;

	private String currentTaxLabel;

	private String currentTreeName;

	private boolean currentTreeStarred;

	private StringBuffer currentTreeRoot = new StringBuffer();

	private StringBuffer currentTreeContent = new StringBuffer();

	/**
	 * Delegates to NexusBlockParser.Abstract.
	 * 
	 * @param blockListener
	 *            the listener to send parse events to.
	 */
	public TreesBlockParser(TreesBlockListener blockListener) {
		super(blockListener);
	}

	public void resetStatus() {
		this.expectingTranslate = true;
		this.expectingTaxLabel = false;
		this.expectingTaxName = false;
		this.expectingTree = true;
		this.expectingTreeName = false;
		this.expectingTreeEquals = false;
		this.expectingTreeRooted = false;
		this.expectingTreeContent = false;
		this.currentTaxLabel = null;
		this.currentTreeStarred = false;
		this.currentTreeName = null;
		this.currentTreeRoot.setLength(0);
		this.currentTreeContent.setLength(0);
	}

	public void parseToken(String token) throws ParseException {
		if (token.trim().length() == 0)
			return;
		else if (this.expectingTranslate && "TRANSLATE".equalsIgnoreCase(token)) {
			this.expectingTranslate = false;
			this.expectingTaxLabel = true;
			this.expectingTree = false;
		} else if (this.expectingTaxLabel) {
			this.currentTaxLabel = token; // In case it includes spaces.
			this.expectingTaxLabel = false;
			this.expectingTaxName = true;
		} else if (this.expectingTaxName) {
			final boolean endsWithComma = token.endsWith(",");
			final String taxName = endsWithComma ? token.substring(0, token
					.length() - 1) : token;
			((TreesBlockListener) this.getBlockListener()).addTranslation(
					this.currentTaxLabel, taxName);
			this.expectingTaxName = false;
			if (!endsWithComma)
				this.expectingTree = true;
			else 
				this.expectingTaxLabel = true;
		} else if (this.expectingTree && "TREE".equalsIgnoreCase(token)) {
			this.expectingTree = false;
			this.expectingTreeName = true;
		} else if (this.expectingTreeName) {
			if ("*".equals(token))
				this.currentTreeStarred = true;
			else {
				this.expectingTreeName = false;
				if (token.indexOf("=") >= 0) {
					final String parts[] = token.split("=");
					this.currentTreeName = parts[0];
					if (parts.length > 1) {
						this.currentTreeContent.append(token);
						this.expectingTreeRooted = false;
					}
					this.expectingTreeContent = true;
				} else {
					this.currentTreeName = token;
					this.expectingTreeEquals = true;
				}
			}
		} else if (this.expectingTreeEquals && token.startsWith("=")) {
			this.expectingTreeEquals = false;
			final String parts[] = token.split("=");
			if (parts.length > 1) {
				this.currentTreeContent.append(parts[1]);
				this.expectingTreeRooted = false;
			} else
				this.expectingTreeRooted = true;
			this.expectingTreeContent = true;
		} else if (this.expectingTreeContent) {
			this.currentTreeContent.append(token); // In case it includes
			// spaces.
			this.expectingTreeRooted = false;
		} else
			throw new ParseException("Found unexpected token " + token
					+ " in TREES block");
	}

	public void endTokenGroup() {
		if (this.expectingTreeContent) {
			final NewickTreeString tree = new NewickTreeString();
			tree.setRootType(this.currentTreeRoot.toString());
			tree.setTreeString(this.currentTreeContent.toString());
			tree.setStarred(this.currentTreeStarred);
			((TreesBlockListener) this.getBlockListener()).addTree(
					this.currentTreeName, tree);
			this.currentTreeContent.setLength(0);
			this.currentTreeName = null;
			this.currentTreeRoot.setLength(0);
			this.currentTreeStarred = false;
			this.expectingTreeRooted = false;
			this.expectingTreeContent = false;
			this.expectingTree = true;
		} else
			super.endTokenGroup();
	}

}
