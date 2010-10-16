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

interface ParentContext extends CodeContext {
    public void promoteOutstandingReference(OutstandingReference ref);
    public void writeShortAt(int pos, int i);
    public int getOffset();

    public int resolveLocalNoCreate(LocalVariable lv);
    public void setMaxLocals(int m);
    public int getUsedLocals();

    public void addExceptionTableEntry(ExceptionMemento em);
}
