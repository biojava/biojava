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


/**
 * Builds a Nexus block from listening to events.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface NexusBlockBuilder extends NexusBlockListener {

	/**
	 * Obtain the constructed block.
	 * 
	 * @return the constructed block.
	 */
	public NexusBlock getNexusBlock();

	/**
	 * This abstract version knows how to build and add comments.
	 */
	public abstract class Abstract implements NexusBlockBuilder {

		private String blockName;

		private NexusBlock block;

		private NexusComment comment;

		/**
		 * Obtains the name of this block.
		 */
		protected String getBlockName() {
			return this.blockName;
		}

		public void beginComment() {
			if (this.comment != null)
				this.comment.openSubComment();
			else
				this.comment = new NexusComment();
		}

		public void commentText(String comment) {
			this.comment.addCommentText(comment);
		}

		public void endComment() {
			if (this.comment != null && this.comment.hasOpenSubComment())
				this.comment.closeSubComment();
			else {
				this.addComment(this.comment);
				this.comment = null;
			}
		}

		/**
		 * Tell the builder to add the given comment at the current location.
		 * 
		 * @param comment
		 *            the comment to add.
		 * @throws ParseException
		 *             if the comment was invalid.
		 */
		protected abstract void addComment(NexusComment comment);

		public void startBlock(String blockName) {
			this.blockName = blockName;
			this.block = this.startBlockObject();
		}

		/**
		 * Tell the builder to start a new block object.
		 */
		protected abstract NexusBlock startBlockObject();

		public NexusBlock getNexusBlock() {
			return this.block;
		}
	}
}
