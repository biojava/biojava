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
 * A CodeClass implementation that is used to generate new classes.
 *
 * <p>
 * When creating classes, instantiate one of these, add fields and methods.
 * Associate CodeGenerator instances with methods. Then, use
 * GeneratedClassLoader to make a new class.
 * </p>
 *
 * @author Matthew Pocock
 */
public class GeneratedCodeClass implements CodeClass {
  private String name;
  private CodeClass superClass;
  private List interfaces;
  private int modifiers;
  private Map methods;
  private Map fields;
  private String sourceFile;
  private boolean deprecated;

  {
    methods = new HashMap();
    fields = new HashMap();
    sourceFile = null;
  }

  public GeneratedCodeClass(
          String name,
          Class superClass,
          Class[] interfaces,
          int modifiers
          ) throws CodeException
  {
    this.name = name;
    this.modifiers = modifiers;
    this.superClass = IntrospectedCodeClass.forClass(superClass);
    this.interfaces = new ArrayList(Arrays.asList(interfaces));
    for (Iterator i = this.interfaces.iterator(); i.hasNext();) {
      Class clazz = (Class) i.next();
      if (!clazz.isInterface()) {
        throw new CodeException(
                "Attempted to create class implemneting non-interface " + clazz
        );
      }
    }
  }

  public GeneratedCodeClass(String name,
                            CodeClass superClass,
                            CodeClass[] interfaces,
                            int modifiers)
          throws CodeException
  {
    this.name = name;
    this.modifiers = modifiers;
    this.superClass = superClass;
    this.interfaces = new ArrayList(Arrays.asList(interfaces));
    for (Iterator i = this.interfaces.iterator(); i.hasNext();) {
      Object obj = i.next();
      if (!(obj instanceof CodeClass)) {
        throw new CodeException(
                "Interface list must contain CodeClass instances"
        );
      }
    }
  }

  /**
   * Set the source file associated with this code class.
   *
   * <p>
   * The source file appears in debugging output and stack traces. Use this
   * method to set the source file that this generated class will clame to be
   * from. You can use non-file names e.g. uri:myGenerator:proxy/foo
   * </p>
   *
   * <p>
   * To un-set the source file, use null.
   * </p>
   *
   * @param sourceFile  the source file for this class
   */
  public void setSourceFile(String sourceFile) {
    this.sourceFile = sourceFile;
  }

  /**
   * Get the source file associated with this code class.
   *
   * <p>
   * Null indicates that no source file is set.
   * </p>
   *
   * @return the source file for this code class
   */
  public String getSourceFile() {
    return sourceFile;
  }

  /**
   * Set the deprecation flag.
   *
   * <p>
   * If deprecated is true, the class will be flagged as deprecated.
   * </p>
   *
   * @param deprecated  the new value of the deprecation
   */
  public void setDeprecated(boolean deprecated) {
    this.deprecated = deprecated;
  }

  /**
   * Get the deprecation flag.
   *
   * @return  wether or not this class is deprecated
   */
  public boolean isDeprecated() {
    return deprecated;
  }

  public List getInterfaces() {
    return Collections.unmodifiableList(interfaces);
  }

  public Set getMethods() {
    return methods.keySet();
  }

  public Set getMethodsByName(String name) {
    Set all = getMethods();
    Set some = new HashSet();
    for (Iterator i = all.iterator(); i.hasNext();) {
      CodeMethod m = (CodeMethod) i.next();
      if (m.getName().equals(name)) {
        some.add(m);
      }
    }
    return some;
  }

  public CodeMethod getConstructor(CodeClass[] args)
          throws NoSuchMethodException
  {
    return getMethod("<init>", args);
  }

  public CodeMethod getMethod(String name, CodeClass[] args)
          throws NoSuchMethodException
  {
    Set poss = getMethodsByName(name);
    METHOD_LOOP:
     for (Iterator i = poss.iterator(); i.hasNext();) {
       CodeMethod meth = (CodeMethod) i.next();
       if (meth.numParameters() != args.length) {
         continue METHOD_LOOP;
       }
       for (int j = 0; j < args.length; j++) {
         if (!meth.getParameterType(j).equals(args[j])) {
           continue METHOD_LOOP;
         }
       }
       return meth;
     }

    StringBuffer methodSig = new StringBuffer(
            "Could not find method " + getName() + "." + name + "("
    );
    if (args.length > 0) {
      methodSig.append(args[0].getName());
    }
    for (int i = 1; i < args.length; i++) {
      methodSig.append(",");
      methodSig.append(args[i].getName());
    }
    methodSig.append(")");
    throw new NoSuchMethodException(methodSig.toString());
  }

  public Set getFields() {
    return fields.keySet();
  }

  public CodeClass getSuperClass() {
    return superClass;
  }

  public CodeField getFieldByName(String name)
          throws NoSuchFieldException
  {
    CodeField f = (CodeField) fields.get(name);
    if (f == null) {
      throw new NoSuchFieldException("No field for " + name + " in class " + getName());
    }
    return f;
  }

  public String getName() {
    return name;
  }

  public String getJName() {
    String name = getName();
    StringBuffer sb = new StringBuffer();
    for (int i = 0; i < name.length(); ++i) {
      char c = name.charAt(i);
      if (c == '.')
        sb.append('/');
      else
        sb.append(c);
    }
    return sb.toString();
  }

  public int getModifiers() {
    return modifiers;
  }

  public String getDescriptor() {
    String name = getName();
    StringBuffer sb = new StringBuffer();
    sb.append('L');
    for (int i = 0; i < name.length(); ++i) {
      char c = name.charAt(i);
      if (c == '.')
        sb.append('/');
      else
        sb.append(c);
    }
    sb.append(';');
    return sb.toString();
  }

  /**
   * Create a new method.
   *
   * <p>
   * This defines the shape of a method that will be generated. Use
   * {@link #setCodeGenerator} to associate code with the method.
   * </p>
   *
   * <p>
   * The argNames will become the names of local variables for each argument.
   * </p>
   *
   * @param name      the method name
   * @param type      the return type
   * @param args      arguments taken
   * @param argNames  names of the arguments
   * @param mods      access modifiers
   * @return          a new GeneratedCodeMethod
   * @throws CodeException if the method could not be created
   */
  public GeneratedCodeMethod createMethod(
          String name,
          CodeClass type,
          CodeClass[] args,
          String[] argNames,
          int mods
          )
          throws CodeException
  {
    GeneratedCodeMethod cm = new GeneratedCodeMethod(this, name, type, args, argNames, mods);
    if (methods.containsKey(cm)) {
      throw new CodeException("Attempt to create multiple methods with same signatures.");
    }

    methods.put(cm, null);
    return cm;
  }

  /**
   * Create a new method.
   *
   * <p>
   * This defines the shape of a method that will be generated. Use
   * {@link #setCodeGenerator} to associate code with the method.
   * </p>
   *
   * @param name      the method name
   * @param type      the return type
   * @param args      arguments taken
   * @param mods      access modifiers
   * @return          a new GeneratedCodeMethod
   * @throws CodeException if the method could not be created
   */
  public GeneratedCodeMethod createMethod(
          String name,
          CodeClass type,
          CodeClass[] args,
          int mods
          )
          throws CodeException
  {
    return createMethod(name, type, args, new String[0], mods);
  }

  public CodeField createField(String name, CodeClass clazz, int mods)
          throws CodeException
  {
    if (fields.containsKey(name)) {
      throw new CodeException("Attempt to create multiple fields named " + name);
    }

    CodeField cf = new CodeField(this, name, clazz, mods);
    fields.put(name, cf);
    return cf;
  }

  public void setCodeGenerator(CodeMethod method, CodeGenerator cg)
          throws CodeException
  {
    if (!methods.containsKey(method)) {
      throw new CodeException("Class doesn't provide method " + method.getName());
    }

    methods.put(method, cg);
  }

  public void createCode(OutputStream os)
          throws IOException, CodeException
  {
    DataOutputStream dos = new DataOutputStream(os);

    // Write classfile header

    dos.writeInt((int) (0xcafebabe));    // Magic
    dos.writeShort(3);                  // Minor  version
    dos.writeShort(45);                   // Major version (check!)

    ConstantPool cp = new ConstantPool();

    // The rest of the classfile gets written to a buffer, accumulating a constant pool along the way

    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    DataOutputStream bdos = new DataOutputStream(baos);

    bdos.writeShort(modifiers);
    bdos.writeShort(cp.resolveClass(this));         // this-class ID
    bdos.writeShort(cp.resolveClass(superClass));   // super-class ID
    bdos.writeShort(interfaces.size());             // number_of_interfaces
    for (Iterator i = interfaces.iterator(); i.hasNext();) {
      bdos.writeShort(cp.resolveClass((CodeClass) i.next())); // interface ID
    }

    // Write the fields

    bdos.writeShort(fields.size());
    for (Iterator i = fields.values().iterator(); i.hasNext();) {
      CodeField cf = (CodeField) i.next();
      bdos.writeShort(cf.getModifiers());
      bdos.writeShort(cp.resolveUtf8(cf.getName()));
      bdos.writeShort(cp.resolveUtf8(cf.getType().getDescriptor()));
      bdos.writeShort(0); // No attributes right now
    }

    // Write the methods (wahey!)

    Set methSet = methods.entrySet();
    bdos.writeShort(methSet.size());
    for (Iterator i = methSet.iterator(); i.hasNext();) {
      Map.Entry me = (Map.Entry) i.next();
      GeneratedCodeMethod cm = (GeneratedCodeMethod) me.getKey();
      CodeGenerator cg = (CodeGenerator) me.getValue();

      bdos.writeShort(cm.getModifiers());                   // access_flags
      bdos.writeShort(cp.resolveUtf8(cm.getName()));        // name_index
      bdos.writeShort(cp.resolveUtf8(cm.getDescriptor()));  // descriptor_index

      // Actually generate the code
      MethodRootContext ctx = new MethodRootContext(this, cm, cp);
      ctx.open();

      LocalVariable thisP = cm.getThis();
      if (thisP != null) {
        // Non-static method
        ctx.resolveLocal(thisP);
      }
      for (int parm = 0; parm < cm.numParameters(); ++parm) {
        ctx.resolveLocal(cm.getVariable(parm));
      }

      cg.writeCode(ctx);
      ctx.close();

      Set thrownExceptions = cm.getThrownExceptions();

      // number of method attirbutes
      int numMethAttrs = 1; // we always have code

      // do we have exceptions?
      if(!thrownExceptions.isEmpty()) {
        numMethAttrs++;
      }

      bdos.writeShort(numMethAttrs);                        // attributes_count

      // start attribute_info for method

      // Code attribute
      List exceptionTable = ctx.getExceptionTable();

      bdos.writeShort(cp.resolveUtf8("Code"));
      bdos.writeInt(12 + ctx.getOffset() + exceptionTable.size() * 8);
      bdos.writeShort(cg.stackDepth());
      bdos.writeShort(ctx.getMaxLocals());
      bdos.writeInt(ctx.getOffset());
      ctx.writeTo(bdos);
      bdos.writeShort(exceptionTable.size());
      for (Iterator ei = exceptionTable.iterator(); ei.hasNext();) {
        ExceptionMemento em = (ExceptionMemento) ei.next();
        if (!(em.isFullyResolved()))
          throw new CodeException("Exception table entry refers to unresolved label");
        bdos.writeShort(em.startHandled.getOffset());
        bdos.writeShort(em.endHandled.getOffset());
        bdos.writeShort(em.handler.getOffset());
        if (em.eClass != null)
          bdos.writeShort(cp.resolveClass(em.eClass));
        else
          bdos.writeShort(0); // For `finally'
      }
      bdos.writeShort(0); // Code has no sub-attributes

      // Exceptions attribute
      if (thrownExceptions.size() > 0) {
        bdos.writeShort(cp.resolveUtf8("Exceptions"));  // attribute_name_index
        bdos.writeInt(2 + thrownExceptions.size() * 2); // attribute_length
        bdos.writeShort(thrownExceptions.size());       // number_of_exceptions
        for (Iterator tei = thrownExceptions.iterator(); tei.hasNext();) {
          CodeClass exClass = (CodeClass) tei.next();
          bdos.writeShort(cp.resolveClass(exClass));    // exception class
        }
      }
    }

    // class-wide attributes
    //
    // currently, these are SourceFile and Deprecated only
    int classAttributes = 0;

    if(sourceFile != null) {
      classAttributes++;
    }
    if(deprecated) {
      classAttributes++;
    }

    bdos.writeShort(classAttributes);  // attributes_count

    // write the source file attribute
    if(sourceFile != null) {
      bdos.writeShort(cp.resolveUtf8("SourceFile"));  // attribute_name_index
      bdos.writeInt(2);                               // attribute_length
      bdos.writeShort(cp.resolveUtf8(sourceFile));    // sourcefile_index
    }

    // write the deprecate attribute
    if(isDeprecated()) {
      bdos.writeShort(cp.resolveUtf8("Deprecated"));  // attribute_name_index
      bdos.writeInt(0);                               // attribute_length
    }

    // All constants will now have been resolved, so we can finally write the cpool

    dos.writeShort(cp.constantPoolSize());
    cp.writeConstantPool(dos);

    // Append the rest of the classfile to the stream

    baos.writeTo(dos);
  }

  public boolean isPrimitive() {
    return false;
  }

  public boolean isArray() {
    return false;
  }
}
