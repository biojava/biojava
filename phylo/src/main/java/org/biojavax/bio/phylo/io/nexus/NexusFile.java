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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Represents Nexus files.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class NexusFile {

	private List objects = new ArrayList();

	/**
	 * Appends an object to the end of the file.
	 * 
	 * @param object
	 *            the NexusObject to append.
	 */
	public void addObject(final NexusObject object) {
		this.objects.add(object);
	}

	/**
	 * Inserts an object at the given position.
	 * 
	 * @param object
	 *            the NexusObject to insert.
	 * @param pos
	 *            the position (0-indexed) to add it at.
	 */
	public void insertObject(final NexusObject object, final int pos) {
		this.objects.add(pos, object);
	}

	/**
	 * Removes an object from the file.
	 * 
	 * @param object
	 *            the NexusObject to remove.
	 */
	public void removeObject(final NexusObject object) {
		this.objects.remove(object);
	}

	/**
	 * Checks to see if we contain an object.
	 * 
	 * @param object
	 *            the NexusObject to check.
	 * @return <tt>true</tt> if we contain it.
	 */
	public boolean containsObject(final NexusObject object) {
		return this.objects.contains(object);
	}

	/**
	 * Iterate over all objects in the file in order.
	 * 
	 * @return an iterator of NexusObjects.
	 */
	public Iterator objectIterator() {
		return this.objects.iterator();
	}

	/**
	 * Iterate over all comments in the file in order.
	 * 
	 * @return an iterator of NexusComments.
	 */
	public Iterator commentIterator() {
		final List comments = new ArrayList();
		for (final Iterator i = this.objectIterator(); i.hasNext();) {
			final NexusObject obj = (NexusObject) i.next();
			if (obj instanceof NexusComment)
				comments.add(obj);
		}
		return comments.iterator();
	}

	/**
	 * Iterate over all blocks in the file in order.
	 * 
	 * @return an iterator of NexusBlocks.
	 */
	public Iterator blockIterator() {
		final List blocks = new ArrayList();
		for (final Iterator i = this.objectIterator(); i.hasNext();) {
			final NexusObject obj = (NexusObject) i.next();
			if (obj instanceof NexusBlock)
				blocks.add(obj);
		}
		return blocks.iterator();
	}
}
