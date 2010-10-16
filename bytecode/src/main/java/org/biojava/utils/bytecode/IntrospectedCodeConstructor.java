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
package org.biojava.utils.bytecode;

import java.lang.reflect.*;

/**
 * Wrap up details about a constructor in a Java class file
 *
 * @author Matthew Pocock
 */

class IntrospectedCodeConstructor implements CodeMethod {
    private final Constructor _constructor;

    public IntrospectedCodeConstructor(Constructor m) {
	_constructor = m;
    }

    public String getName() {
	return "<init>"; // all constructors are called <init>
    }

    public String getFullName() {
	return _constructor.getDeclaringClass().getName() + "." + getName();
    }

    public CodeClass getContainingClass() {
	return IntrospectedCodeClass.forClass(_constructor.getDeclaringClass());
    }

    public String getDescriptor() {
	StringBuffer sb = new StringBuffer();
	sb.append('(');
	for (int i = 0; i < numParameters(); ++i) {
	    CodeClass cc = getParameterType(i);
	    sb.append(cc.getDescriptor());
	}
	sb.append(')');
	sb.append(getReturnType().getDescriptor());
	return sb.toString();
    }

    public int getModifiers() {
	return _constructor.getModifiers();
    }

    public CodeClass getReturnType() {
	return CodeUtils.TYPE_VOID;
    }

    public int numParameters() {
	return _constructor.getParameterTypes().length;
    }

    public CodeClass getParameterType(int pos) {
	return IntrospectedCodeClass.forClass(_constructor.getParameterTypes()[pos]);
    }
}
