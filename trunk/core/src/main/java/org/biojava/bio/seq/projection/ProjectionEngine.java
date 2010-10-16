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

package org.biojava.bio.seq.projection;

import java.io.InputStream;
import java.lang.reflect.Constructor;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.Feature;
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
import org.biojava.utils.bytecode.InstructionVector;
import org.biojava.utils.bytecode.IntrospectedCodeClass;
import org.biojava.utils.bytecode.LocalVariable;

/**
 * Factory for proxy objects which project BioJava features
 * into alternate coordinate systems.  This class binds
 * together a feature and a <code>ProjectionContext</code>
 * object, and returns the resulting projected feature.
 * New feature-projection wrapper classes are generated
 * automatically as they are required.
 *
 * Don't use this class directly. This class contains deep voodoo code. Run
 * away while you still can.
 *
 * You may find that for some bizaare reason you need to manually project a feature
 * through a ProjectionContext. You should be using @link
 * ProjectionContext.projectFeature(). If this is not practical for some reason,
 * use this class. However, this probably indicates that you are doing something
 * mad.
 *
 * Projected feature classes will be named as
 * org.biojava.bio.seq.projection.{$ctxt}.{$feat} where $ctxt and $feat are the
 * full class names of the context and feature respectively with each "."
 * character replaced with "_".
 *
 * Delegate into this from your ProjectionContext to do the dirty work of
 * actually projecting features from the underlying data. This factory
 * generate a new class that is unique to the combination of your
 * projection context and the Feature interface that the projected feature
 * implements. Then, it will instantiate that class with your context and the
 * underlying feature. The returned feature will currently be an implementation
 * of the package-private class org.biojava.bio.seq.projection.ProjectedFeature,
 * but this is implementation detail and should not be relied upon.
 *
 * @author Thomas Down
 * @since 1.2
 */

public class ProjectionEngine {
  /**
   * The standard projection engine object.
   */

  public static final ProjectionEngine DEFAULT;

  static {
    DEFAULT = new ProjectionEngine();
  }

  private ProjectionEngine() {
    super();
  }

  private final Map _projectionClasses;
  private final PEClassLoader loader;
  private final Instantiator instantiator;

  {
    loader = new PEClassLoader(ClassTools.getClassLoader(ProjectionEngine.class));
    _projectionClasses = new HashMap();
    try {
      instantiator = new InstantiatorImpl();
    } catch (Exception ex) {
      throw new AssertionFailure(
              "Assertion failure: can't initialize projection system", ex);
    }
  }

  private Class searchForClass(Class origClass, Class ctxtClass) {
    Map pc = (Map) _projectionClasses.get(ctxtClass);

    if(pc == null) {
      return null;
    }

    return (Class) pc.get(origClass);
  }

  private void registerClass(Class origClass,
                             Class ctxtClass,
                             Class projClass)
  {
    Map pc = (Map) _projectionClasses.get(ctxtClass);

    if(pc == null) {
      pc = new HashMap();
      pc.put(null, ProjectedFeature.class);
      _projectionClasses.put(ctxtClass, pc);
    }

    pc.put(origClass, projClass);
  }

  /**
   * Return a projection of Feature <code>f</code> into the system
   * defined by a given ProjectionContext.  The returned object
   * will implement the same Feature interface (sub-interface of
   * <code>Feature</code> as the underlying feature, and will also
   * implement the <code>Projection</code> interface.
   */

  public Feature projectFeature(Feature f, ProjectionContext ctx) {
    //System.err.println("Searching with : " + f + " : " + ctx);
    Class featureClass = f.getClass();
    Class[] fcInterfaces = featureClass.getInterfaces();
    Class featureInterface = Feature.class;
    for (int i = fcInterfaces.length - 1; i >= 0; --i) {
      if (Feature.class.isAssignableFrom(fcInterfaces[i])) {
        featureInterface = fcInterfaces[i];
        break;
      }
    }

    Class projectionClass = getFeatureProjectionClass(featureInterface,
                                                      ctx.getClass());
    //System.err.println("Got projectionClass: " + projectionClass + " for " + featureInterface + " : " + ctx.getClass());
    Class[] sig = new Class[2];
    sig[0] = featureInterface;
    sig[1] = ProjectionContext.class;
    try {
      Constructor ct = projectionClass.getConstructor(sig);
      Object[] args = new Object[2];
      args[0] = f;
      args[1] = ctx;
      return (Feature) instantiator.newInstance(ct, args);
    } catch (Exception ex) {
      throw new AssertionFailure(
              "Assertion failed: Couldn't instantiate proxy " +
              projectionClass.getName(),
              ex);
    }
  }

  private synchronized Class getFeatureProjectionClass(Class face, Class ctxtClass) {
    Class projection = searchForClass(face, ctxtClass);
    if (projection == null) {
      try {
        Class baseClass = ProjectedFeature.class;

        String faceName = face.getName()
                .replaceAll("\\.", "_").replaceAll("\\$", "__");
        String ctxtName = ctxtClass.getName()
                .replaceAll("\\$", "_");

        CodeClass baseClassC = IntrospectedCodeClass.forClass(baseClass);
        CodeClass faceClassC = IntrospectedCodeClass.forClass(face);

        GeneratedCodeClass pclass = new GeneratedCodeClass(
                ctxtName + "__" + faceName,
                baseClassC,
                new CodeClass[]{faceClassC},
                CodeUtils.ACC_PUBLIC | CodeUtils.ACC_SUPER
        );
        //System.err.println("Generating proxy class: " + pclass.getName());

        CodeClass[] baseInitArgsList = new CodeClass[] {
          IntrospectedCodeClass.forClass(Feature.class),
          IntrospectedCodeClass.forClass(ProjectionContext.class)
        };

        CodeClass voidC
                = IntrospectedCodeClass.forClass(Void.TYPE);
        CodeClass projectionContextC
                = IntrospectedCodeClass.forClass(ProjectionContext.class);
        CodeClass ourContextC
                = IntrospectedCodeClass.forClass(ctxtClass);

        CodeMethod m_ourBase_init = baseClassC.getConstructor(baseInitArgsList);

        CodeMethod m_ourBase_getViewedFeature = baseClassC.getMethod(
                "getViewedFeature",
                CodeUtils.EMPTY_LIST);

        CodeMethod m_ourBase_getProjectionContext = baseClassC.getMethod(
                "getProjectionContext",
                CodeUtils.EMPTY_LIST);

        GeneratedCodeMethod init = pclass.createMethod(
                "<init>",
                voidC,
                new CodeClass[]{
                  faceClassC,
                  projectionContextC
                },
                CodeUtils.ACC_PUBLIC);

        InstructionVector initIV = new InstructionVector();
        initIV.add(ByteCode.make_aload(init.getThis()));
        initIV.add(ByteCode.make_aload(init.getVariable(0)));
        initIV.add(ByteCode.make_aload(init.getVariable(1)));
        initIV.add(ByteCode.make_invokespecial(m_ourBase_init));
        initIV.add(ByteCode.make_return());
        pclass.setCodeGenerator(init, initIV);

        METHOD_MAKER:
        for (Iterator methIt = faceClassC.getMethods().iterator(); methIt.hasNext();) {
          CodeMethod faceMethod = (CodeMethod) methIt.next();
          Set baseMethods = baseClassC.getMethodsByName(faceMethod.getName());

          if (baseClassC.getMethodsByName(faceMethod.getName()).size() > 0) {
            for(Iterator i = baseMethods.iterator(); i.hasNext(); ) {
              CodeMethod meth = (CodeMethod) i.next();
              if( (meth.getModifiers() & CodeUtils.ACC_ABSTRACT) == 0) {
                //System.err.println("Skipping defined method: " + faceMethod.getName());
                continue METHOD_MAKER;
              }
            }
          }

          //System.err.println("Looking at method: " + faceMethod.getName());
          if (faceMethod.getName().startsWith("get") &&
                  faceMethod.numParameters() == 0) {
            //System.err.println("Getter: " + faceMethod.getName());
            String propName = faceMethod.getName().substring("get".length());

            // we will make a proxy
            GeneratedCodeMethod getterMethod = pclass.createMethod(
                    faceMethod.getName(),
                    faceMethod.getReturnType(),
                    CodeUtils.EMPTY_LIST,
                    CodeUtils.ACC_PUBLIC);

            // search for project*
            //System.err.println("Searching for methods: project" + propName + " : " + ourContextC.getMethodsByName("project" + propName));
            CodeMethod projMeth = null;
            for(Iterator i = ourContextC.getMethodsByName(
                    "project" + propName).iterator();
                i.hasNext();) {
              CodeMethod cm = (CodeMethod) i.next();
              //System.err.println("Evaluating context method: " + cm.getName() + " return: " + cm.getReturnType().getName() + " params: " + cm.numParameters());
              if(cm.numParameters() == 1 &&
                 cm.getReturnType().equals(cm.getParameterType(0)) &&
                 cm.getReturnType().equals(faceMethod.getReturnType())
              ) {
                //System.err.println("A match");
                projMeth = cm;
                break;
              }
            }

            if(projMeth == null) {
              // direct proxy
              //
              // equivalent to:
              //  Foo ProjFeat.getFoo() {
              //    return getViewedFeature().getFoo();
              //  }

              InstructionVector proxyIV = new InstructionVector();
              proxyIV.add(ByteCode.make_aload(getterMethod.getThis()));
              proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getViewedFeature));
              proxyIV.add(ByteCode.make_invokeinterface(faceMethod));
              proxyIV.add(ByteCode.make_return(getterMethod));
              pclass.setCodeGenerator(getterMethod, proxyIV);
              //System.err.println("Direct proxy: " + getterMethod);
            } else {
              // context proxy
              //
              // equivalent to:
              //   Foo projFeat.getFoo() {
              //      return getContext().projectFoo(
              //         getViewedFeature().getFoo());
              //  }

              InstructionVector proxyIV = new InstructionVector();
              proxyIV.add(ByteCode.make_aload(getterMethod.getThis()));
              proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getProjectionContext));
              proxyIV.add(ByteCode.make_checkcast(ourContextC));
              proxyIV.add(ByteCode.make_aload(getterMethod.getThis()));
              proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getViewedFeature));
              proxyIV.add(ByteCode.make_invokeinterface(faceMethod));
              proxyIV.add(ByteCode.make_invokevirtual(projMeth));
              proxyIV.add(ByteCode.make_return(getterMethod));
              pclass.setCodeGenerator(getterMethod, proxyIV);
              //System.err.println("TargetContext proxy: " + getterMethod + " " + projMeth);
            }
          }

          if (faceMethod.getName().startsWith("set") &&
                  faceMethod.numParameters() == 0) {
            // System.err.println("Setter: " + faceMethod.getName());
            String propName = faceMethod.getName().substring("set".length());

            // we will make a proxy
            GeneratedCodeMethod setterMethod = pclass.createMethod(
                    faceMethod.getName(),
                    voidC,
                    new CodeClass[] { faceMethod.getParameterType(0) },
                    CodeUtils.ACC_PUBLIC);

            // search for revert*
            CodeMethod revertMeth = null;
            for(Iterator i = ourContextC.getMethodsByName(
                    "revert" + propName).iterator();
                i.hasNext();) {
              CodeMethod cm = (CodeMethod) i.next();
              if(cm.numParameters() == 1 &&
                 cm.getReturnType().equals(cm.getParameterType(0)) &&
                 cm.getReturnType().equals(faceMethod.getReturnType())
              ) {
                revertMeth = cm;
                break;
              }

              if(revertMeth == null) {
                // direct proxy
                //
                // equivalent to:
                //  void ProjFeat.setFoo(Foo foo) {
                //    getViewedFeature().setFoo(foo);
                //  }

                InstructionVector proxyIV = new InstructionVector();
                proxyIV.add(ByteCode.make_aload(setterMethod.getThis()));
                proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getViewedFeature));
                proxyIV.add(ByteCode.make_lload(setterMethod.getVariable(0)));
                proxyIV.add(ByteCode.make_invokeinterface(faceMethod));
                pclass.setCodeGenerator(setterMethod, proxyIV);
                //System.err.println("Direct proxy: " + setterMethod);
              } else {
                // context proxy
                //
                // equivalent to:
                //  void ProjFeat.setFoo(Foo foo) {
                //    getViewedFeature().setFoo(
                //      getProjectionContext().revertFoo(foo));
                //  }
                InstructionVector proxyIV = new InstructionVector();
                proxyIV.add(ByteCode.make_aload(setterMethod.getThis()));
                proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getViewedFeature));
                proxyIV.add(ByteCode.make_aload(setterMethod.getThis()));
                proxyIV.add(ByteCode.make_invokevirtual(m_ourBase_getProjectionContext));
                proxyIV.add(ByteCode.make_checkcast(ourContextC) );
                proxyIV.add(ByteCode.make_lload(setterMethod.getVariable(0)));
                proxyIV.add(ByteCode.make_invokevirtual(revertMeth));
                proxyIV.add(ByteCode.make_invokeinterface(faceMethod));
                pclass.setCodeGenerator(setterMethod, proxyIV);
                //System.err.println("TargetContext proxy: " + setterMethod + " " + revertMeth);
              }
            }
          }
        }

        //System.err.println("Creating class: " + pclass);
        projection = loader.defineClass(pclass);
        registerClass(face, ctxtClass, projection);
      } catch (CodeException ex) {
        throw new BioError(ex);
      } catch (NoSuchMethodException nsme) {
        throw new BioError(nsme);
      }

    }
    return projection;
  }

  /**
   * Revert a template so that it can be used on the original feature-space.
   *
   * <p>
   * This will use the revertFoo methods defined in the context to revert
   * all of the transformed properties.
   * </p>
   *
   * @param templ the template to revert
   * @param ctxt  the context defining the reversion
   * @return  a new template of the same type as templ, with reverted fields
   */
  public Feature.Template revertTemplate(Feature.Template templ, ProjectionContext ctxt) {
    Class templateClass = templ.getClass();
    Class projClass = getTemplateProjectorClass(templateClass, ctxt.getClass());
    TemplateProjector proj = null;
    try {
      proj = (TemplateProjector) projClass.newInstance();
    } catch (InstantiationException ie) {
      throw new AssertionFailure(
              "Assertion failed: Couldn't instantiate template projector" +
              projClass.getName(),
              ie);
    } catch (IllegalAccessException iae) {
      throw new AssertionFailure(
              "Assertion Failed: Couldn't instantiate template projector" +
              projClass.getName(),
              iae);
    }

    try {
      return proj.revertTemplate(ctxt, templ);
    } catch (CloneNotSupportedException cnse) {
      throw new AssertionFailure(cnse);
    }
  }

  private Class getTemplateProjectorClass(Class templateClass, Class ctxtClass) {
    Class projection = searchForClass(templateClass, ctxtClass);
    if(projection == null) {
      try {
        String tpltName = templateClass.getName()
                .replaceAll("\\.", "_").replaceAll("\\$", "__");
        String ctxtName = ctxtClass.getName()
                .replaceAll("\\$", "_");

        CodeClass c_baseClass = CodeUtils.TYPE_OBJECT;
        CodeClass c_tpClass = IntrospectedCodeClass.forClass(TemplateProjector.class);

        GeneratedCodeClass pclass = new GeneratedCodeClass(
                ctxtName + "__" + tpltName,
                c_baseClass,
                new CodeClass[] { c_tpClass },
                CodeUtils.ACC_PUBLIC | CodeUtils.ACC_SUPER);

        // we need a dumb no-args constructor
        CodeMethod m_Object_init
                = CodeUtils.TYPE_OBJECT.getConstructor(CodeUtils.EMPTY_LIST);
        GeneratedCodeMethod init = pclass.createMethod("<init>",
                                                       CodeUtils.TYPE_VOID,
                                                       CodeUtils.EMPTY_LIST,
                                                       CodeUtils.ACC_PUBLIC);
        InstructionVector initIV = new InstructionVector();
        initIV.add(ByteCode.make_aload(init.getThis()));
        initIV.add(ByteCode.make_invokespecial(m_Object_init));
        initIV.add(ByteCode.make_return());
        pclass.setCodeGenerator(init, initIV);

        // we need to implement revertTemplate
        CodeClass c_FeatureTemplate = IntrospectedCodeClass.forClass(Feature.Template.class);
        CodeClass c_ProjectionContext = IntrospectedCodeClass.forClass(ProjectionContext.class);
        GeneratedCodeMethod revT = pclass.createMethod("revertTemplate",
                                                       c_FeatureTemplate,
                                                       new CodeClass[] {
                                                         c_ProjectionContext,
                                                         c_FeatureTemplate },
                                                       new String[] {
                                                         "ctxt", "origTplt" },
                                                       CodeUtils.ACC_PUBLIC);
        revT.addThrownException(IntrospectedCodeClass.forClass(CloneNotSupportedException.class));
        InstructionVector revTIV = new InstructionVector();

        // firstly, revertTempalte must cast the arg to the correct type and
        // then clone it and store this away in a local variable
        CodeMethod m_FeatureTemplate_clone = c_FeatureTemplate.getMethod(
                "clone",
                CodeUtils.EMPTY_LIST);
        LocalVariable lv_ctxt = revT.getVariable("ctxt");
        LocalVariable lv_origTplt = revT.getVariable("origTplt");
        LocalVariable lv_ourTplt = new LocalVariable(c_FeatureTemplate, "ourTemplate");
        revTIV.add(ByteCode.make_aload(lv_origTplt));
        revTIV.add(ByteCode.make_invokevirtual(m_FeatureTemplate_clone));
        revTIV.add(ByteCode.make_checkcast(c_FeatureTemplate));
        revTIV.add(ByteCode.make_astore(lv_ourTplt));

        // This method must assign each field that is mapped using revertFoo.
        // The pattern will be for each revertable field, call
        // ((TemplateClass) ourTemplate).foo =
        //    ((TemplateClass) origTemplate).foo
        CodeClass c_ourContext = IntrospectedCodeClass.forClass(ctxtClass);
        CodeClass c_ourTemplate = IntrospectedCodeClass.forClass(templateClass);

        for(Iterator fi = c_ourTemplate.getFields().iterator();
            fi.hasNext(); ) {
          CodeField field = (CodeField) fi.next();
          String fieldName = field.getName();
          String propName = fieldName.substring(0, 1).toUpperCase() +
                  fieldName.substring(1);
          String revName = "revert" + propName;

          CodeMethod revertMeth = null;
          for(Iterator mi = c_ourContext.getMethodsByName(
                  revName).iterator();
              mi.hasNext(); )
          {
            CodeMethod cm = (CodeMethod) mi.next();
            if(cm.numParameters() == 1 &&
                    cm.getReturnType().equals(cm.getParameterType(0)) &&
                    cm.getReturnType().equals(field.getType())) {
              revertMeth = cm;
              break;
            }
          }

          if(revertMeth != null) {
            // we have a revert method. Wire it in using the pattern:
            //
            // ((TemplateClass) ourTemplate).foo =
            //    ((OurCtxt) ctxt).revertFoo(
            //      ((TemplateClass) origTemplate).foo)
            InstructionVector revIV = new InstructionVector();
            revIV.add(ByteCode.make_aload(lv_ourTplt));
            revIV.add(ByteCode.make_checkcast(c_ourTemplate));
            revIV.add(ByteCode.make_aload(lv_ctxt));
            revIV.add(ByteCode.make_checkcast(c_ourContext));
            revIV.add(ByteCode.make_aload(lv_origTplt));
            revIV.add(ByteCode.make_checkcast(c_ourTemplate));
            revIV.add(ByteCode.make_getfield(field));
            revIV.add(ByteCode.make_invokevirtual(revertMeth));
            revIV.add(ByteCode.make_putfield(field));
            revTIV.add(revIV);
          }
        }
        revTIV.add(ByteCode.make_aload(lv_ourTplt));
        revTIV.add(ByteCode.make_areturn());

        pclass.setCodeGenerator(revT, revTIV);

        projection = loader.defineClass(pclass);
        registerClass(templateClass, ctxtClass, projection);
      } catch (CodeException ce) {
        throw new AssertionFailure("Unable to create template projector: ", ce);
      } catch (NoSuchMethodException nsme) {
        throw new AssertionFailure("Unable to create template projector: ", nsme);
      }
    }
    return projection;
  }

  private static class PEClassLoader extends GeneratedClassLoader {
    public PEClassLoader(ClassLoader parent) {
      super(parent);
    }

    public Class loadClassMagick(String name)
            throws ClassNotFoundException {
      try {
        String resName = name.replace('.', '/') + ".class";
        InputStream is = getResourceAsStream(resName);
        if (is == null) {
          throw new ClassNotFoundException(
                  "Could not find class resource: " + resName);
        }

        byte[] buffer = new byte[10000];
        int len = 0;
        int read = 0;
        while (read >= 0) {
          read = is.read(buffer, len, buffer.length - len);
          if (read > 0) {
            len += read;
          }
        }
        is.close();
        Class c = defineClass(name, buffer, 0, len);
        resolveClass(c);

        return c;
      } catch (Exception ex) {
        throw new ClassNotFoundException(
                "Could not load class for " + name + ": " + ex.toString());
      }

    }
  }

  /**
   * Internal helper class.
   */

  public static interface Instantiator {
    public Object newInstance(Constructor c, Object[] args) throws Exception;
  }

  /**
   * Internal helper class.
   */

  static class InstantiatorImpl implements Instantiator {
    public Object newInstance(Constructor c, Object[] args) throws Exception {
      return c.newInstance(args);
    }
  }

  /**
   * This is an interface for things that project feature templates.
   *
   * <p><em>Note:</em> This is implementation guts that has to be public to bind
   * contexts to features. It is not intended for users.</p>
   *
   * @author Matthew Pocock
   * @since 1.4
   */
  public static interface TemplateProjector {
    public Feature.Template revertTemplate(ProjectionContext ctxt,
                                           Feature.Template projTempl)
            throws CloneNotSupportedException;
  }
}
