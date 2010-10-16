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
 * A method that will be generated.
 *
 * <p>
 * These are instantiated by factory methods on {@link GeneratedCodeClass}, and
 * can not be instantiated directly.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public final class GeneratedCodeMethod implements CodeMethod {
  private final String name;
  private List args;
  private List localvars;
  private LocalVariable thisV;
  private int modifiers;
  private CodeClass type;
  private CodeClass container;
  private Set thrownExceptions;
  private Map nameToLocals;

  GeneratedCodeMethod(
    CodeClass container, 
    String name, 
    CodeClass type, 
    CodeClass[] args, 
    String[] names,
    int modifiers
  ) {
    this.container = container;
    this.name = name;
    this.args = new ArrayList(Arrays.asList(args));
    this.modifiers = modifiers;
    this.type = type;
    nameToLocals = new HashMap();
    localvars = new ArrayList();
    for(int i = 0; i < this.args.size(); ++i) {
      if(i < names.length) {
        LocalVariable arg = new LocalVariable(args[i], names[i]);
        localvars.add(arg);
        nameToLocals.put(names[i], arg);
      } else {
        localvars.add(new LocalVariable(args[i]));
      }
    }

    if((modifiers & CodeUtils.ACC_STATIC) == 0) {
      thisV = new LocalVariable(container, "this");
      nameToLocals.put("this", thisV);
    }
  }

  public String getName() {
    return name;
  }

  public String getFullName() {
    return container.getName() + "." + name;
  }

  public CodeClass getContainingClass() {
    return container;
  }

  public String getDescriptor() {
    StringBuffer sb = new StringBuffer();
    sb.append('(');
    for(Iterator i = args.iterator(); i.hasNext(); ) {
      CodeClass cc = (CodeClass) i.next();
      sb.append(cc.getDescriptor());
    }
    sb.append(')');
    sb.append(type.getDescriptor());
    return sb.toString();
  }

  public int getModifiers() {
    return modifiers;
  }

  public CodeClass getReturnType() {
    return type;
  }

  public int numParameters() {
    return args.size();
  }

  public CodeClass getParameterType(int pos) {
    return (CodeClass) args.get(pos);
  }

  /**
   *  Gets the Variable attribute of the GeneratedCodeMethod object.
   *
   * <p>
   * There is one local variable for each of the arguments of the method,
   * indexed from 0.
   * </p>
   *
   * @param  pos  the index of the local variable
   * @return      the local variable
   */
  public LocalVariable getVariable(int pos) {
    return (LocalVariable) localvars.get(pos);
  }
  
  /**
   * Gets the Variable attribute of the GenerateCodeMethod object by name.
   *
   * <P>
   * All methods have a variable under the string "this". If it was constructed
   * with a String [] naming the args, the locals for each local can be
   * retrieved by name.
   * </p>
   *
   * @param argName a String naming the local
   * @return        the LocalVariable for that argName
   * @throws        NoSuchElementException if there is no local with that name
   */
  public LocalVariable getVariable(String argName)
  throws NoSuchElementException {
    LocalVariable lv = (LocalVariable) nameToLocals.get(argName);
    if(lv == null) {
      throw new NoSuchElementException(
        "Can't find local for argName " + argName
      );
    }
    return lv;
  }

  /**
   *  Gets the This attribute of the GeneratedCodeMethod object 
   *
   * @return    The This value 
   */
  public LocalVariable getThis() {
    return thisV;
  }

  /**
   *  Gets the ThrownExceptions attribute of the GeneratedCodeMethod object 
   *
   * @return    The ThrownExceptions value 
   */
  public Set getThrownExceptions() {
    return Collections.unmodifiableSet(thrownExceptions);
  }

  /**
   *  Adds a feature to the ThrownException attribute of the GeneratedCodeMethod object 
   *
   * @param  cc  The feature to be added to the ThrownException attribute 
   */
  public void addThrownException(CodeClass cc) {
    thrownExceptions.add(cc);
  }
  
  {
    thrownExceptions = new HashSet();
  }
}
