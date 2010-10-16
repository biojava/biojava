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
import java.io.*;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 */
class MethodRootContext implements CodeContext, ParentContext {
  private final CodeClass codeClass;
  private final CodeMethod codeMethod;
  private final ConstantPool cpool;
  
  private byte[] sink;
  private int offset;
  
  private List outstandingRefs;
  private Map markedLabels;
  
  private Map localVariables;
  
  private int usedLocals;
  private int maxLocals;
  private Map resolvedParametrics;

  private List exceptionTable;
  
  {
    outstandingRefs = new ArrayList();
    markedLabels = new HashMap();
    localVariables = new HashMap();
    exceptionTable = new ArrayList();
    resolvedParametrics = new HashMap();
  }
  
  MethodRootContext(CodeClass cc, CodeMethod cm, ConstantPool cp) {
    this.codeClass = cc;
    this.codeMethod = cm;
    this.cpool = cp;
    
    sink = new byte[65536];
    offset = 0;
    
    // Should allocate local variables from method;
    
    usedLocals = maxLocals = 0;
  }
  
  public CodeClass getCodeClass() {
    return codeClass;
  }
  
  public CodeMethod getCodeMethod() {
    return codeMethod;
  }
  
  public ConstantPool getConstants() {
    return cpool;
  }
  
  public void writeByte(byte b) {
    sink[offset++] = b;
  }
  
  public void writeShort(int i) {
    sink[offset++] = (byte) ((i>>8) & 0xff);
    sink[offset++] = (byte) ((i) & 0xff);
  }
  
  public void writeShortAt(int pos, int i) {
    sink[pos] = (byte) ((i>>8) & 0xff);
    sink[pos+1] = (byte) ((i) & 0xff);
  }
  
  public void writeLabel(Label l) {
    outstandingRefs.add(new BranchFixup(l, getOffset(), this));
    writeShort(0);
  }
  
  public void markLabel(Label l) throws CodeException {
    if (markedLabels.containsKey(l))
      throw new CodeException("Attempt to duplicate marked label");
    markedLabels.put(l, new Integer(getOffset()));
  }
  
  public int resolveLocal(LocalVariable lv) {
    // System.out.println("MethodRootContext.resolveLocal(" + lv + ")");
    Integer slot = (Integer) localVariables.get(lv);
    if (slot != null) {
      // System.out.println("resolved to " + slot);
      return slot.intValue();
    }
    int newLocal = usedLocals;
    usedLocals += lv.needSlots();
    setMaxLocals(usedLocals); 
    localVariables.put(lv, new Integer(newLocal));
    // System.out.println("created at" + newLocal);
    return newLocal;
  }
  
  public int resolveLocalNoCreate(LocalVariable lv) {
    // System.out.println("MethodRootContext.resolveLocalNoCreate(" + lv + ")");
    Integer slot = (Integer) localVariables.get(lv);
    if (slot != null) {
      // System.out.println("at " + slot);
      return slot.intValue();
    } else {
      // System.out.println("not here");
      return -1;
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
      throw new CodeException("Can not resolve type: " + type);
    }
    
    return cc;
  }
  
  public int getOffset() {
    return offset;
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
        throw new CodeException("Reference to label " + or.getLabel() + " still outstanding at top level");
      }
    }
  }
  
  public CodeContext subContext() {
    return new ChildContext(this);
  }
  
  public void promoteOutstandingReference(OutstandingReference or) {
    outstandingRefs.add(or);
  }
  
  public void writeTo(OutputStream os) throws IOException {
    os.write(sink, 0, offset);
  }
  
  public void setMaxLocals(int newMax) {
    if (newMax > maxLocals) 
      maxLocals = newMax;
  }
  
  public int getMaxLocals() {
    return maxLocals;
  }
  
  public int getUsedLocals() {
    return usedLocals;
  }
  
  public void addExceptionTableEntry(
    Label startHandled, 
    Label endHandled,
    CodeClass eClass, 
    Label handler
  ) {
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
    exceptionTable.add(em);
  }
  
  public List getExceptionTable() {
    return Collections.unmodifiableList(exceptionTable);
  }
}
