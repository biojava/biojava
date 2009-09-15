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

import java.io.IOException;
import java.io.Writer;

/**
 * Represents a Nexus block.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface NexusBlock extends NexusObject {

	/**
	 * Get the block name.
	 * 
	 * @return the block name.
	 */
	public String getBlockName();

	public abstract class Abstract implements NexusBlock {
		private String blockName;

		/**
		 * Construct a block with a given name.
		 * 
		 * @param blockName
		 *            the name to give the block.
		 */
		public Abstract(final String blockName) {
			this.blockName = blockName;
		}

		public String getBlockName() {
			return this.blockName;
		}

		public void writeObject(final Writer writer) throws IOException {
			writer.write("BEGIN " + this.blockName + ";");
			writer.write(NexusFileFormat.NEW_LINE);
			this.writeBlockContents(writer);
			writer.write("END;");
			writer.write(NexusFileFormat.NEW_LINE);
		}

		/**
		 * Writes a token and correctly substitutes all symbols in it.
		 * 
		 * @param writer
		 *            the writer to write to.
		 * @param token
		 *            the token to write.
		 * @throws IOException
		 *             if writing failed.
		 */
		protected void writeToken(final Writer writer, String token)
				throws IOException {
			token = token.replaceAll("'", "''");
			token = token.replaceAll("_", "'_'");
			if (token.trim().length() > 0)
				token = token.replaceAll(" ", "_");
			if (token.indexOf('[') >= 0 || token.indexOf(']') >= 0)
				token = "'" + token + "'";
			writer.write(token);
		}

		/**
		 * Implement this to write out block contents, not including the BEGIN
		 * and END tags.
		 * 
		 * @param writer
		 *            the writer to write to.
		 * @throws IOException
		 *             if writing failed.
		 */
		protected abstract void writeBlockContents(Writer writer)
				throws IOException;
	}
}
