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

import java.util.*;

/**
 * CodeContext implementation which provides a lightweight subContext of any
 * ParentContext.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class ChildContext implements CodeContext, ParentContext {
  ParentContext ourParent;
  
  private List outstandingRefs;
  private Map markedLabels;
  
  private Map localVariables;
  private int usedLocals;
  
  private Map resolvedParametrics;
  
  {
    outstandingRefs = new ArrayList();
    markedLabels = new HashMap();
    localVariables = new HashMap();
    resolvedParametrics = new HashMap();
  }

  ChildContext(ParentContext p) {
    this.ourParent = p;
    usedLocals = ourParent.getUsedLocals();
  }

  public CodeClass getCodeClass() {
    return ourParent.getCodeClass();
  }
  
  public CodeMethod getCodeMethod() {
    return ourParent.getCodeMethod();
  }
  
  public ConstantPool getConstants() {
    return ourParent.getConstants();
  }
  
  public void writeByte(byte b) throws CodeException {
    ourParent.writeByte(b);
  }
  
  public void writeShort(int i) throws CodeException {
    ourParent.writeShort(i);
  }
  
  public void writeShortAt(int pos, int i) {
    ourParent.writeShortAt(pos, i);
  }
  
  public void markLabel(Label l) throws CodeException {
    if (markedLabels.containsKey(l))
      throw new CodeException("Attempt to duplicate marked label");
    markedLabels.put(l, new Integer(getOffset()));
  }
  
  public void promoteOutstandingReference(OutstandingReference or) {
    outstandingRefs.add(or);
  }
  
  public void writeLabel(Label l) throws CodeException {
    outstandingRefs.add(new BranchFixup(l, getOffset(), this));
    writeShort(0);
  }
  
  public int getOffset() {
    return ourParent.getOffset();
  }
  
  public void open() {
  }
  
  public void close() throws CodeException {
    for (ListIterator li = outstandingRefs.listIterator(); li.hasNext(); ) {
      OutstandingReference or = (OutstandingReference) li.next();
      Integer off = (Integer) markedLabels.get(or.getLabel());
      if (off != null) {
        or.resolve(off.intValue());
      } else {
        ourParent.promoteOutstandingReference(or); // See if parent can resolve this.
      }
    }
  }
  
  public CodeContext subContext() {
    return new ChildContext(this);
  }
  
  public int resolveLocal(LocalVariable lv) {
    // System.out.println("ChildContext.resolveLocal(" + lv + ")");
    Integer slot = (Integer) localVariables.get(lv);
    if (slot != null) {
      // System.out.println("resolved to " + slot);
      return slot.intValue();
    }
    
    // System.out.println("Trying to resolve slot in parent " + ourParent);
    int locSlot = ourParent.resolveLocalNoCreate(lv);
    if (locSlot >= 0) {
      // System.out.println("Parent had local " + lv + " at " + locSlot);
    } else {
      // Need to create the variable;
      locSlot = usedLocals;
      usedLocals += lv.needSlots();
      setMaxLocals(usedLocals);
      // System.out.println("Generated new slot for local " + lv + " at " + locSlot);
    }
    
    // We'll add the slot to our map, even if it's just a copy from the parent.
    
    localVariables.put(lv, new Integer(locSlot));
    // System.out.println("ChildContext: Resolved local variable " + lv + " to slot " + locSlot);
    return locSlot;
  }
  
  public int resolveLocalNoCreate(LocalVariable lv) {
    Integer slot = (Integer) localVariables.get(lv);
    if (slot != null) {
      return slot.intValue();
    } else {
      return ourParent.resolveLocalNoCreate(lv); 
    }
  }
  
  public void registerParametricType(
    ParametricType type,
    CodeClass concreteType
  ) throws CodeException {
    if(resolvedParametrics.containsKey(type)) {
      throw new CodeException("Failed to regiter parametric type " + type +
        ". Attempted to register for " + concreteType +
        " but it is already registered for " + resolvedParametrics.get(type) );
    }
    
    if(!type.canAccept(concreteType)) {
      throw new CodeException(
        "Parametric type is not compattible with concrete type: " +
        type + " : " + concreteType );
    }

    resolvedParametrics.put(type, concreteType);
  }
  
  public CodeClass resolveParametricType(ParametricType type)
  throws CodeException {
    CodeClass cc = (CodeClass) resolvedParametrics.get(type);
    
    if(cc == null) {
      return ourParent.resolveParametricType(type);
    }
    
    return cc;
  }
  
  public int getUsedLocals() {
    return usedLocals;
  }
  
  public void setMaxLocals(int m) {
    ourParent.setMaxLocals(m);
  }
  
  public void addExceptionTableEntry(Label startHandled, 
  Label endHandled,
  CodeClass eClass, 
  Label handler)
  {
    SimpleReference rStartHandled = new SimpleReference(startHandled);
    SimpleReference rEndHandled = new SimpleReference(endHandled);
    SimpleReference rHandler = new SimpleReference(handler);
    outstandingRefs.add(rStartHandled);
    outstandingRefs.add(rEndHandled);
    outstandingRefs.add(rHandler);
    
    addExceptionTableEntry(new ExceptionMemento(rStartHandled,
    rEndHandled,
    eClass,
    rHandler));
  }
  
  public void addExceptionTableEntry(ExceptionMemento em) {
    ourParent.addExceptionTableEntry(em);
  }
}
