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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Represents a Nexus comment, possibly nested.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class NexusComment implements NexusObject {

	private List components = new ArrayList();

	private NexusComment nestedComment = null;

	public void openSubComment() {
		if (this.hasOpenSubComment())
			this.nestedComment.openSubComment();
		else {
			this.nestedComment = new NexusComment();
			this.components.add(this.nestedComment);
		}
	}

	public boolean hasOpenSubComment() {
		return this.nestedComment != null;
	}

	public void closeSubComment() {
		if (this.hasOpenSubComment() && this.nestedComment.hasOpenSubComment())
			this.nestedComment.closeSubComment();
		else
			this.nestedComment = null;
	}

	public void addCommentText(final String text) {
		if (this.hasOpenSubComment())
			this.nestedComment.addCommentText(text);
		else
			this.components.add(text);
	}

	/**
	 * This iterator iterates over all parts of the comment. Each item returned
	 * is either a String or a NexusComment.
	 * 
	 * @return an iterator over the comment components.
	 */
	public Iterator commentIterator() {
		return this.components.iterator();
	}

	public void writeObject(final Writer writer) throws IOException {
		writer.write('[');
		for (final Iterator i = this.components.iterator(); i.hasNext();) {
			final Object obj = (Object) i.next();
			if (obj instanceof NexusComment)
				((NexusComment) obj).writeObject(writer);
			else {
				String text = (String) obj;
				text = text.replaceAll("'", "''");
				if (text.indexOf('[') >= 0 || text.indexOf(']') >= 0)
					text = "'" + text + "'";
				writer.write(text);
			}
		}
		writer.write(']');
	}

}
