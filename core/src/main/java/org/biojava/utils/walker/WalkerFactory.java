/**
 * Factory for generating walkers that are customised to particuar feature
 * visitors.
 *
 * @author Matthew Pocock
 */
package org.biojava.utils.walker;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ClassTools;
import org.biojava.utils.bytecode.ByteCode;
import org.biojava.utils.bytecode.CodeClass;
import org.biojava.utils.bytecode.CodeException;
import org.biojava.utils.bytecode.CodeField;
import org.biojava.utils.bytecode.CodeMethod;
import org.biojava.utils.bytecode.CodeUtils;
import org.biojava.utils.bytecode.GeneratedClassLoader;
import org.biojava.utils.bytecode.GeneratedCodeClass;
import org.biojava.utils.bytecode.GeneratedCodeMethod;
import org.biojava.utils.bytecode.IfExpression;
import org.biojava.utils.bytecode.InstructionVector;
import org.biojava.utils.bytecode.IntrospectedCodeClass;
import org.biojava.utils.bytecode.LocalVariable;

public class WalkerFactory {
  private static Map factories = new HashMap();

  /**
   * Make a WalkerFactory that handles a Visitor for
   * a class of type typeClazz.
   *
   * @param typeClazz  the Class this factory will walk over
   */
  public synchronized static WalkerFactory getInstance(Class typeClazz)
  {
    WalkerFactory instance;
    String typeName = typeClazz.getName();

    if ((instance = (WalkerFactory) factories.get(typeName))  == null) {
      instance = new WalkerFactory(typeClazz);
      factories.put(typeName, instance);
    }
    return instance;
  }

  public synchronized static WalkerFactory getInstance() {
    return getInstance(FeatureFilter.class);
  }

  private final Map walkers;
  private final GeneratedClassLoader classLoader;
  private final List typesWithParents;
  private final Class typeClazz;

  private WalkerFactory(Class typeClazz) {
    walkers = new HashMap();
    classLoader = new GeneratedClassLoader(ClassTools.getClassLoader(this));
    typesWithParents = new ArrayList();
    this.typeClazz = typeClazz;
  }

  public Class getTypeClass()
  {
    return typeClazz;
  }

  /**
   * Register a type as being a 'container' class.
   * Container classes will be scanned for methods for retrieving child
   * instances that can be walked to.
   *
   * 
   * You should never need to call this. The library authors should take care
   * of this for you.
   *
   * 
   * Register 'structural' classes here - those with children.
   *
   * @param type  the Class of the type with children
   */
  public synchronized void addTypeWithParent(Class type) {
    if (!typesWithParents.contains(type)) {
        typesWithParents.add(type);
    }
  }

  /**
   * Get a Walker that is customosed to a particular visitor.
   *
   * @param visitor  the Visitor this walker will scan with
   * @return  a Walker bound to this visitor
   * @throws BioException if the walker could not be built
   */
  public synchronized Walker getWalker(Visitor visitor)
  throws BioException {
    Class walkerClass = (Class) walkers.get(visitor.getClass());

    if(walkerClass == null) {
      walkers.put(visitor.getClass(),
                  walkerClass = generateWalker(visitor.getClass()));
    }

    try {
      return (Walker) walkerClass.newInstance();
    } catch (InstantiationException ie) {
      throw new AssertionFailure("Could not instantiate walker for class: " +
                                 walkerClass,
                                 ie);
    } catch (IllegalAccessException iae) {
      throw new AssertionFailure("Could not instantiate walker for class: " +
                                 walkerClass,
                                 iae);
    }
  }

  private Class generateWalker(Class visitorClass)
  throws BioException {
    try {
      String vcn = visitorClass.getName().replaceAll("\\$", "_");
      String walkerClassName = vcn + "_walker";

      CodeClass c_Visitor = IntrospectedCodeClass.forClass(Visitor.class);
      CodeClass c_ourVisitor = IntrospectedCodeClass.forClass(visitorClass);
      CodeClass c_Walker = IntrospectedCodeClass.forClass(Walker.class);
      CodeClass c_WalkerBase = IntrospectedCodeClass.forClass(Object.class);
      CodeClass c_Object = IntrospectedCodeClass.forClass(Object.class);

      // make a class
      GeneratedCodeClass walkerClass = new GeneratedCodeClass(
              walkerClassName,
              c_WalkerBase,
              new CodeClass[]{c_Walker},
              CodeUtils.ACC_PUBLIC | CodeUtils.ACC_SUPER);

      // get all visitor methods
      Method[] methods = visitorClass.getMethods();
      List visitorMeths = new ArrayList();
      Class retClass = null;

      // make a set of all handlers and get the return type used
      // barf if the return type is not consistent
      for(int mi = 0; mi < methods.length; mi++) {
        Method method = methods[mi];

        Class ret = method.getReturnType();
        Class[] args = method.getParameterTypes();

        // is this one of our classes?
        if(args.length > 0) {
          Class arg0 = args[0];

          // If arg0 is inner class, strip off outer-class name to make our name
          String methName = method.getName();
          String arg0Name = arg0.getName();
          int doli = arg0Name.lastIndexOf('$');
          if(doli >= 0) {
            arg0Name = arg0Name.substring(doli+1);
          }
          doli = arg0Name.lastIndexOf('.');
          if(doli >= 0) {
            arg0Name = arg0Name.substring(doli+1);
          }

          // drop the leading captial
          arg0Name = arg0Name.substring(0,1).toLowerCase() +
                  arg0Name.substring(1);

          // we have a naming match?
          if(arg0Name.equals(methName)) {
            // check argument 0 is of the type we are visiting
            if(typeClazz.isAssignableFrom(arg0)) {
              // we have a live one.
              // check that the return type is good
              if(retClass == null) {
                retClass = ret;
              } else {
                if(retClass != ret) {
                  throw new BioException(
                          "Return type of all methods must agree. " +
                          "We were expecting: " + retClass.getName() +
                          " but found: " + ret.getName());
                }
              }

              // if there are other args, make sure they match the return type
              for(int ai = 1; ai < args.length; ai++) {
                Class argI = args[ai];
                if(argI != retClass) {
                  throw new BioException(
                          "The first argument to a handler method must be a " +
                          typeClazz.toString() +
                          "All subsequent arguments must match the return type.  In: " +
                          method);
                }
              }

              // OK - this looks like a good handler - add it to our list
              visitorMeths.add(method);
            }
          }
        }
      }

      // if retclass was never set, it's safe to make it void
      if(retClass == null) {
        retClass = Void.TYPE;
      }

      // sort by type - most derived first
      Collections.sort(visitorMeths,  new Comparator() {
        public int compare(Object o1, Object o2) {
          Method m1 = (Method) o1; Class c1 = m1.getParameterTypes()[0];
          Method m2 = (Method) o2; Class c2 = m2.getParameterTypes()[0];
          if(c1.isAssignableFrom(c2)) return +1;
          if(c2.isAssignableFrom(c1)) return -1;
          return 0;
        }
      });

      // now let's implement our dispatcher
      CodeClass c_retClass = IntrospectedCodeClass.forClass(retClass);
      CodeClass[] walkParams = new CodeClass[]{ c_Object, c_Visitor};

      GeneratedCodeMethod doWalk = walkerClass.createMethod(
              "doWalk",
              c_retClass,
              walkParams,
              new String[]{"target", "visitor"},
              CodeUtils.ACC_PUBLIC);
      InstructionVector walkIV = new InstructionVector();
      LocalVariable lv_target = doWalk.getVariable("target");
      LocalVariable lv_visitor = doWalk.getVariable("visitor");
      LocalVariable lv_visitor2 = new LocalVariable(c_ourVisitor, "visitor2");
      walkIV.add(ByteCode.make_aload(lv_visitor));
      walkIV.add(ByteCode.make_checkcast(c_ourVisitor));
      walkIV.add(ByteCode.make_astore(lv_visitor2));

      // local variables for the return values of wrapped invocatins
      List wrappedLVs = new ArrayList();



      // firstly, we should call And, Or, Not, etc., wrapped targets
      //
      // These are all listed in typesWithParents.
      InstructionVector wrapperIV = new InstructionVector();
      for(Iterator fwpi = typesWithParents.iterator(); fwpi.hasNext(); ) {
        InstructionVector wfiv = new InstructionVector();

        // find the methods that get the wrapped
        Class filtClass = (Class) fwpi.next();
        CodeClass c_ourType = IntrospectedCodeClass.forClass(filtClass);

        Method[] filtMeth = filtClass.getMethods();
        int lvi = 0;

        for(int mi = 0; mi < filtMeth.length; mi++) {
          Method m = filtMeth[mi];

          // no args, returns a typeClazz
          if(m.getParameterTypes().length == 0 &&
                  typeClazz.isAssignableFrom(m.getReturnType()))
          {
            CodeMethod m_getChild = IntrospectedCodeClass.forMethod(m);
            LocalVariable lv = null;

            if (c_retClass != CodeUtils.TYPE_VOID) {
              if (lvi < wrappedLVs.size()) {
                lv = (LocalVariable) wrappedLVs.get(lvi);
              } else {
                lv = new LocalVariable(c_retClass);
                wrappedLVs.add(lv);
              }
              lvi++;
            }

            // res_i = this.walk(
            //       ((c_ourType) target).m_getChild(),
            //       visitor );
            wfiv.add(ByteCode.make_aload(doWalk.getThis()));
            wfiv.add(ByteCode.make_aload(lv_target));
            wfiv.add(ByteCode.make_checkcast(c_ourType));
            wfiv.add(ByteCode.make_invokevirtual(m_getChild));
            wfiv.add(ByteCode.make_aload(lv_visitor));
            wfiv.add(ByteCode.make_invokevirtual(doWalk));

            if(c_retClass != CodeUtils.TYPE_VOID) {
              wfiv.add(ByteCode.make_astore(lv));
            }
          }
        }

        // if (target instanceof ourType) {
        //   do the above block
        // }
        wrapperIV.add(ByteCode.make_aload(lv_target));
        wrapperIV.add(ByteCode.make_instanceof(c_ourType));
        wrapperIV.add(new IfExpression(ByteCode.op_ifne,
                                    wfiv,
                                    CodeUtils.DO_NOTHING));
      }

      for(Iterator lvi = wrappedLVs.iterator(); lvi.hasNext(); ) {
        LocalVariable lv = (LocalVariable) lvi.next();
        walkIV.add(ByteCode.make_aconst_null());
        walkIV.add(ByteCode.make_astore(lv));
      }
      walkIV.add(wrapperIV);

      // the big if/else/if/else stack goes here - switching on the
      // type for each method using instanceof
      //
      // if(filter instanceof Filt1) {
      //    viewer2.filt1( ((Filt1) filter) );
      //    return;
      // }

      for(Iterator mi = visitorMeths.iterator(); mi.hasNext(); ) {
        Method meth = (Method) mi.next();

        // the viewer method is invoked as:
        //
        //   viewer2.foo( (Foo) filter, ...);
        //
        // if the return value is void, we just return
        // if it is not void, we return the Bar instance it returns

        CodeMethod ourMeth = IntrospectedCodeClass.forMethod(meth);
        CodeClass c_thisFiltType = ourMeth.getParameterType(0);

        InstructionVector bodyIV = new InstructionVector();
        bodyIV.add(ByteCode.make_aload(lv_visitor2));
        bodyIV.add(ByteCode.make_aload(lv_target));
        bodyIV.add(ByteCode.make_checkcast(c_thisFiltType));
        for(int ai = 1; ai < ourMeth.numParameters(); ai++) {
          bodyIV.add(ByteCode.make_aload((LocalVariable) wrappedLVs.get(ai-1)));
        }
        bodyIV.add(ByteCode.make_invokevirtual(ourMeth));
        bodyIV.add(ByteCode.make_return(doWalk));

        // wrap this in an if(filt instanceof Foo)
        walkIV.add(ByteCode.make_aload(lv_target));
        walkIV.add(ByteCode.make_instanceof(c_thisFiltType));
        walkIV.add(new IfExpression(ByteCode.op_ifne,
                                    bodyIV,
                                    CodeUtils.DO_NOTHING));
      }

      // return void if we are void,
      // return null if we are meant to return something but no handler was used
      if(c_retClass == CodeUtils.TYPE_VOID) {
        walkIV.add(ByteCode.make_return());
      } else {
        walkIV.add(ByteCode.make_aconst_null());
        walkIV.add(ByteCode.make_areturn());
      }

      walkerClass.setCodeGenerator(doWalk, walkIV);

      // Wire doWalk to walk and if necisary create field
      //
      CodeField f_value = null;
      if(c_retClass != CodeUtils.TYPE_VOID) {
        f_value = walkerClass.createField("value",
                                          CodeUtils.TYPE_OBJECT,
                                          CodeUtils.ACC_PRIVATE);
      }

      { // protect us from leakey locals
        GeneratedCodeMethod walkImpl = walkerClass.createMethod(
                "walk",
                CodeUtils.TYPE_VOID,
                walkParams,
                CodeUtils.ACC_PUBLIC);
        InstructionVector wiIV = new InstructionVector();
        if (c_retClass != CodeUtils.TYPE_VOID) {
          wiIV.add(ByteCode.make_aload(walkImpl.getThis()));
        }
        wiIV.add(ByteCode.make_aload(walkImpl.getThis()));
        wiIV.add(ByteCode.make_aload(walkImpl.getVariable(0)));
        wiIV.add(ByteCode.make_aload(walkImpl.getVariable(1)));
        wiIV.add(ByteCode.make_invokevirtual(doWalk));
        if (c_retClass != CodeUtils.TYPE_VOID) {
          wiIV.add(ByteCode.make_putfield(f_value));
        }
        wiIV.add(ByteCode.make_return());

        walkerClass.setCodeGenerator(walkImpl, wiIV);
      }

      // generate the getValue() method
      { // protect us from leakey locals
        GeneratedCodeMethod getValue = walkerClass.createMethod(
                "getValue",
                CodeUtils.TYPE_OBJECT,
                CodeUtils.EMPTY_LIST,
                CodeUtils.ACC_PUBLIC);
        InstructionVector gvIV = new InstructionVector();
        if (c_retClass == CodeUtils.TYPE_VOID) {
          gvIV.add(ByteCode.make_aconst_null());
        } else {
          gvIV.add(ByteCode.make_aload(getValue.getThis()));
          gvIV.add(ByteCode.make_getfield(f_value));
        }
        gvIV.add(ByteCode.make_areturn());
        walkerClass.setCodeGenerator(getValue, gvIV);
      }

      // constructor - no args, forward to SUPER, intialize field if needed
      { // protect us from leaky locals
        CodeMethod m_WalkerBase_init = c_WalkerBase.getConstructor(CodeUtils.EMPTY_LIST);
        GeneratedCodeMethod init = walkerClass.createMethod("<init>",
                                                            CodeUtils.TYPE_VOID,
                                                            CodeUtils.EMPTY_LIST,
                                                            CodeUtils.ACC_PUBLIC);
        InstructionVector initIV = new InstructionVector();
        initIV.add(ByteCode.make_aload(init.getThis()));
        initIV.add(ByteCode.make_invokespecial(m_WalkerBase_init));
        if (c_retClass != CodeUtils.TYPE_VOID) {
          initIV.add(ByteCode.make_aload(init.getThis()));
          initIV.add(ByteCode.make_aconst_null());
          initIV.add(ByteCode.make_putfield(f_value));
        }
        initIV.add(ByteCode.make_return());
        walkerClass.setCodeGenerator(init, initIV);
      }

      return classLoader.defineClass(walkerClass);
    } catch (CodeException ce) {
      throw new AssertionFailure("Unable to generate walker code", ce);
    } catch (NoSuchMethodException nsme) {
      throw new AssertionFailure("Unable to generate walker code", nsme);
    }
  }
}

