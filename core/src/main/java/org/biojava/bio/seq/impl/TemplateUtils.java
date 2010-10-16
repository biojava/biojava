package org.biojava.bio.seq.impl;

import java.lang.reflect.Field;
import java.lang.reflect.Method;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.utils.AssertionFailure;

/**
 * Common things you may want to do with feature templates.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public final class TemplateUtils {
  // no instances of this class
  private TemplateUtils() {}

  /**
   * This attempts to divine the 'best' template class for a feature and return
   * a new instance readly for pupulating.
   *
   * <p>
   * This will (hopefully) be the most derived feature interface implemented by
   * a feature class. This code assumes that feature interfaces are singly
   * inherited. Of course, with the current template system, it is a fairly safe
   * assumption.
   * </p>
   *
   *
   */
  public static Feature.Template instantiateTemplate(Feature feat)
  throws BioException {
    Feature.Template tmpl = searchForTemplate(feat.getClass());
    if(tmpl == null) {
      throw new AssertionFailure(
              "Could not find template for feature class: " + feat.getClass());
    }
    return tmpl;
  }

  private static Feature.Template searchForTemplate(Class fClass)
  throws BioException {
    if(fClass.isInterface()) {
      if(Feature.class.isAssignableFrom(fClass)) {
        return instantiateTemplate(fClass);
      }
    }

    Class[] interfaces = fClass.getInterfaces();
    for(int i = 0; i < interfaces.length; i++) {
      Feature.Template tmpl = searchForTemplate(interfaces[i]);
      if(tmpl != null) {
        return tmpl;
      }
    }

    Class parent = fClass.getSuperclass();
    if(parent != null) {
      return searchForTemplate(parent);
    }

    return null;
  }

  private static Feature.Template instantiateTemplate(Class fClass)
  throws BioException{
    // let's assume this has a *.Template class & instantiate that
    Class[] declClasses = fClass.getDeclaredClasses();
    for(int i = 0; i < declClasses.length; i++) {
      if(declClasses[i].getName().equals(fClass.getName() + "$" + "Template")) {
        try {
          return (Feature.Template) declClasses[i].newInstance();
        } catch (IllegalAccessException iae) {
          throw new BioException(
                  "Expected template no-args constructor to be accessible",
                  iae);
        } catch (InstantiationException ie) {
          throw new BioException(
                  "Failed to execute no-args constructor",
                  ie);
        }
      }
    }

    return null;
  }

  /**
   * This attempts to populate the fields of a template using
   * the publically accessible information in a feature. It is simple
   * to call populate() within Feature.makeTemplate() to ensure all the
   * slots get filled.
   *
   * @param templ the Feature.Template instance to populate
   * @param feat  the Feature to read info from
   */
  public static void populate(Feature.Template templ, Feature feat)
  throws BioException {
    Field[] fields = templ.getClass().getFields();
    Method[] methods = feat.getClass().getMethods();

    for(int i = 0; i < fields.length; i++) {
      Field field = fields[i];
      String fName = field.getName();
      String methName =
        "get" +
        fName.substring(0, 1).toUpperCase() +
        fName.substring(1);

      Method method = null;
      for(int j = 0; j < methods.length; j++) {
        Method meth = methods[j];
        if(methods[j].getName().equals(methName)) {
          method = meth;
        }
      }
      if(method == null) {
        throw new BioException("Expecting to find a method named: " + methName);
      }

      try {
        field.set(templ, method.invoke(feat, new Object[] {}));
      } catch (Exception e) {
        throw new BioError("Couldn't access template fields", e);
      }
    }
  }

  public static Feature.Template makeTemplate(Feature feat)
  throws BioException {
    Feature.Template templ = instantiateTemplate(feat);
    populate(templ, feat);
    return templ;
  }
}
