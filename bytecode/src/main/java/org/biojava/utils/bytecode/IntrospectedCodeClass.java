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
import java.lang.reflect.*;

// fix for array types so that their descriptors are not prefixed with L

/**
 * CodeClass instances that represent normal Java Class objects.
 *
 * <p>
 * Instances of IntrospectedCodeClass are generated using the static factory
 * methods named forClass(). These methods ensure that the same
 * IntrospectedCodeClass instance is returned for multiple invocations with
 * the same argument.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class IntrospectedCodeClass implements CodeClass {
  private static Map introspectedClasses;
  private static Map primitiveDescriptors;

  static {
    introspectedClasses = new HashMap();

    primitiveDescriptors = new HashMap();
    primitiveDescriptors.put(Byte.TYPE, "B");
    primitiveDescriptors.put(Character.TYPE, "C");
    primitiveDescriptors.put(Double.TYPE, "D");
    primitiveDescriptors.put(Float.TYPE, "F");
    primitiveDescriptors.put(Integer.TYPE, "I");
    primitiveDescriptors.put(Long.TYPE, "J");
    primitiveDescriptors.put(Short.TYPE, "S");
    primitiveDescriptors.put(Boolean.TYPE, "Z");
    primitiveDescriptors.put(Void.TYPE, "V");

  }

  /**
   * Get the CodeClass for a Java Class.
   *
   * @param c  the Java Class to reflect
   * @return  a CodeClass representing the class
   */
  public static CodeClass forClass(Class c) {
    CodeClass cc = (CodeClass) introspectedClasses.get(c);
    if (cc == null) {
      cc = new IntrospectedCodeClass(c);
      introspectedClasses.put(c, cc);
    }
    return cc;
  }

  /**
   * Get the CodeClass for a Java class name.
   *
   * @param name  the Java class name to reflect
   * @return  a CodeClass representing the class
   */
  public static CodeClass forClass(String name) throws ClassNotFoundException {
    Class c = ClassLoader.getSystemClassLoader().loadClass(name);

    return forClass(c);
  }

  public static CodeMethod forMethod(Method method) {
    return new IntrospectedCodeMethod(method);
  }

  //
  // Instance
  //

  private Class clazz;

  private IntrospectedCodeClass(Class c) {
    this.clazz = c;
  }

  public String getName() {
    return clazz.getName();
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

  public String getDescriptor() {
    if (clazz.isPrimitive()) {
      String desc = (String) primitiveDescriptors.get(clazz);
      if (desc == null) {
        throw new RuntimeException("Unknown primitive type " + clazz.getName() + ", eeek!");
      }
      return desc;
    }

    if (clazz.isArray()) {
      return "[" + IntrospectedCodeClass.forClass(clazz.getComponentType()).getDescriptor();
    } else {
      String name = getName();
      StringBuffer sb = new StringBuffer();
      sb.append('L');
      for (int i = 0; i < name.length(); ++i) {
        char c = name.charAt(i);
        if (c == '.') {
          sb.append('/');
        } else {
          sb.append(c);
        }
      }
      sb.append(';');
      return sb.toString();
    }
  }

  public int getModifiers() {
    return clazz.getModifiers();
  }

  public CodeClass getSuperClass() {
    return IntrospectedCodeClass.forClass(clazz.getSuperclass());
  }

  public List getInterfaces() {
    Class[] interfaces = clazz.getInterfaces();
    return Arrays.asList(interfaces);
  }

  private Set _methods;
  private Map _methsByNameSig;
  private Map _methsByName;

  public Set getMethods() {
    initMethods();
    return _methods;
  }

  private void initMethods() {
    if (_methods == null) {
      Map meths = new HashMap();
      popMeths(this.clazz, meths);
      popIMeths(this.clazz, meths);
      _methods = new HashSet(meths.values());
      _methsByNameSig = new HashMap();
      _methsByName = new HashMap();
      for(Iterator i = _methods.iterator(); i.hasNext(); ) {
        CodeMethod m = (CodeMethod) i.next();
        Set mbn = (Set) _methsByName.get(m.getName());
        if(mbn == null) {
          _methsByName.put(m.getName(), mbn = new HashSet());
        }
        mbn.add(m);
        _methsByNameSig.put(makeNameSig(m), m);
      }
    }
  }

  private void popMeths(Class clazz, Map meths) {
    Method[] methods = clazz.getDeclaredMethods();
    for(int i = 0; i < methods.length; i++) {
      Method m = methods[i];
      ArrayList sigList = new ArrayList();
      sigList.add(m.getName());
      sigList.addAll(Arrays.asList(m.getParameterTypes()));
      if(!meths.containsKey(sigList)) {
        meths.put(sigList, new IntrospectedCodeMethod(m));
      }
    }

    Class sup = clazz.getSuperclass();
    if(sup != null) {
      popMeths(sup, meths);
    }
  }

  private void popIMeths(Class clazz, Map meths) {
    if(clazz.isInterface()) {
      Method[] methods = clazz.getDeclaredMethods();
      for(int i = 0; i < methods.length; i++) {
        Method m = methods[i];
        ArrayList sigList = new ArrayList();
        sigList.add(m.getName());
        sigList.addAll(Arrays.asList(m.getParameterTypes()));
        if(!meths.containsKey(sigList)) {
          meths.put(sigList, new IntrospectedCodeMethod(m));
        }
      }
      Class[] interfaces = clazz.getInterfaces();
      for(int i = 0; i < interfaces.length; i++) {
        popIMeths(interfaces[i], meths);
      }
    }

    Class sup = clazz.getSuperclass();
    if(sup != null) {
      popIMeths(sup, meths);
    }
  }

  private List makeNameSig(CodeMethod m) {
    List res = new ArrayList();
    res.add(m.getName());

    for(int i = 0; i < m.numParameters(); i++) {
      res.add(m.getParameterType(i));
    }

    return res;
  }

  public CodeField getFieldByName(String name)
          throws NoSuchFieldException {
    try {
      Field f = clazz.getField(name);
      return new CodeField(this,
                           name,
                           IntrospectedCodeClass.forClass(f.getType()),
                           f.getModifiers());
    } catch (NoSuchFieldException ex) {
      throw (NoSuchFieldException) new NoSuchFieldException(
              "Can't find field " + name +
              " in class " + getName()
      ).initCause(ex);
    }
  }

  private Set _fields;

  public Set getFields() {
    if(_fields == null) {
      _fields = new HashSet();
      Field[] fields = clazz.getFields();
      for(int fi = 0; fi < fields.length; fi++) {
        Field f = fields[fi];
        _fields.add(new CodeField(this,
                                  f.getName(),
                                  IntrospectedCodeClass.forClass(f.getType()),
                                  f.getModifiers()));
      }

      _fields = Collections.unmodifiableSet(_fields);
    }

    return _fields;
  }

  public Set getMethodsByName(String name) {
    initMethods();

    Set hits = (Set) _methsByName.get(name);
    if(hits == null) {
      return Collections.EMPTY_SET;
    } else {
      return hits;
    }
  }

  public CodeMethod getMethod(String name, CodeClass[] args)
          throws NoSuchMethodException
  {
    initMethods();

    List nameSig = new ArrayList();
    nameSig.add(name);
    for(int i = 0; i < args.length; i++) {
      nameSig.add(args[i]);
    }

    CodeMethod cm = (CodeMethod) _methsByNameSig.get(nameSig);

    if(cm == null) {
      throw new NoSuchMethodException(
              "Could not find method " + getName() +
              "." + name +
              "(" + CodeUtils.classListToString(args) + ")");
    }

    return cm;
  }


  public CodeMethod getConstructor(CodeClass[] args)
          throws NoSuchMethodException {
    try {
      Class[] argsC = new Class[args.length];
      for (int i = 0; i < args.length; i++) {
        argsC[i] = ((IntrospectedCodeClass) args[i]).clazz;
      }
      return new IntrospectedCodeConstructor(clazz.getConstructor(argsC));
    } catch (NoSuchMethodException nsme) {
      throw (NoSuchMethodException) new NoSuchMethodException(
              "Could not find constructor new " + getName() +
              "(" + CodeUtils.classListToString(args) + ")"
      ).initCause(nsme);
    }
  }

  public boolean isPrimitive() {
    return clazz.isPrimitive();
  }

  public boolean isArray() {
    return clazz.isArray();
  }

  public String toString() {
    return this.getClass().getName() + ": " + clazz.getName();
  }
}
