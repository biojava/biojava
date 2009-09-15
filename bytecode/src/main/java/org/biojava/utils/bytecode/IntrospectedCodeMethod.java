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
 * Wrap up details about a method in a Java class file
 *
 * @author Thomas Down
 */

class IntrospectedCodeMethod implements CodeMethod {
    private final Method _method;

    public IntrospectedCodeMethod(Method m) {
	_method = m;
    }

    public String getName() {
	return _method.getName();
    }

    public String getFullName() {
	return _method.getDeclaringClass().getName() + "." + getName();
    }

    public CodeClass getContainingClass() {
	return IntrospectedCodeClass.forClass(_method.getDeclaringClass());
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
	return _method.getModifiers();
    }

    public CodeClass getReturnType() {
	return IntrospectedCodeClass.forClass(_method.getReturnType());
    }

    public int numParameters() {
	return _method.getParameterTypes().length;
    }

    public CodeClass getParameterType(int pos) {
	return IntrospectedCodeClass.forClass(_method.getParameterTypes()[pos]);
    }
  
  public String toString() {
    return super.toString() + " " + getDescriptor();
  }
}
