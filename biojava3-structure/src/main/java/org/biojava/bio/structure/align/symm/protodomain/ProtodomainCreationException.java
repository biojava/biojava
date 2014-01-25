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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

/**
 * An error in the creation of a Protodomain.
 * 
 * @author dmyersturnbull
 */
@SuppressWarnings("serial")
public class ProtodomainCreationException extends Exception {

	private final String protodomainString;
	private final String enclosing;

	private static String makeMsg(String protodomainString, String enclosing, String reason) {
		return "Could not create protodomain " + protodomainString + " (scop ID " + enclosing + "). Reason: " + reason;
	}

	public ProtodomainCreationException(String protodomainString, String enclosing, String reason) {
		super(makeMsg(protodomainString, enclosing, reason));
		this.protodomainString = protodomainString;
		this.enclosing = enclosing;
	}

	public ProtodomainCreationException(String protodomainString, String enclosing, Throwable e) {
		super(makeMsg(protodomainString, enclosing, e.getMessage()), e);
		this.protodomainString = protodomainString;
		this.enclosing = enclosing;
	}

	public ProtodomainCreationException(String protodomainString, String enclosing, Throwable e, String reason) {
		super(makeMsg(protodomainString, enclosing, reason), e);
		this.protodomainString = protodomainString;
		this.enclosing = enclosing;
	}

	public String getProtodomainString() {
		return protodomainString;
	}

	public String getEnclosingName() {
		return enclosing;
	}

}
