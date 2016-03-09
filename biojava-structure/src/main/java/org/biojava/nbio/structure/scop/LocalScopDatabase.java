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
package org.biojava.nbio.structure.scop;

/**
 * Classes which implement ScopDatabase in a way which allows them to serve queries
 * without accessing the internet should implement this interface instead. An
 * initial file download is acceptable.
 *
 * ScopFactory utilizes this distinction to optimize when many queries are expected
 * or remote calls would be otherwise undesireable.
 * @author Spencer Bliven
 *
 */
public interface LocalScopDatabase extends ScopDatabase {}
