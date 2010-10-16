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

package org.biojava.utils.xml;

import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.SmallMap;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * Construct java beans from XML elements
 *
 * @author Thomas Down
 */

public class XMLBeans {
    public final static XMLBeans INSTANCE;

    static {
	INSTANCE = new XMLBeans();
    }

    protected XMLBeans() {
    }

    public Object instantiateBean(Element bel) 
        throws AppException
    {
	return instantiateBean(bel, ClassLoader.getSystemClassLoader(), new HashMap());
    }

    public Object instantiateBean(Element bel, ClassLoader cloader, Map beanRefs) 
        throws AppException
    {
	String cl = bel.getAttribute("jclass");
	if (cl == null)
	    throw new AppException("No jclass attribute");

	Object bean = null;

	try {
	    Class clazz = cloader.loadClass(cl);
	    bean = clazz.newInstance();
	    configureBean(bean, bel, beanRefs);
	    if (bean instanceof Initializable)
		((Initializable) bean).init();   // FIXME
	} catch (ClassNotFoundException ex) {
	    throw new AppException("Couldn't load bean class " + cl);
        } catch (ClassCastException ex) {
	    throw new AppException("Does not implement AppBean: " + cl);
	} catch (InstantiationException ex) {
	    throw new AppException("Couldn't intantiate bean " + cl);
	} catch (IllegalAccessException ex) {
	    throw new AppException("Couldn't access constructor for bean " + cl);
	}

	return bean;
    }

    private void configureBean(Object bean, Element el, Map refs) 
        throws AppException
    {
	Class clazz = bean.getClass();

	Node child = el.getFirstChild();
	while (child != null) {
	    if (child instanceof Element) {
		Element echild = (Element) child;
		String tag = echild.getTagName();
		String name = echild.getAttribute("name");
		Object valueObject = null;
		Class valueType = null;

		if (tag.equals("string")) {
		    valueObject = echild.getAttribute("value");
		    valueType = valueObject.getClass();
		} else if (tag.equals("bean") || tag.equals("child")) {
		    // child supported for backwards compatibility.

		    String ref = echild.getAttribute("ref");
		    Object targ = null;
		    if (! ref.equals("")) {
			targ = refs.get(ref);
                        if(targ == null) {
                          throw new NullPointerException(
                            "Can't find target for: " + ref
                          );
                        }
		    } else {
			targ = instantiateBean(echild);
		    }
		    
		    valueObject = targ;
		    valueType = targ.getClass();
		} else if (tag.equals("int")) {
		    String value = echild.getAttribute("value");
		    try {
			int val = Integer.parseInt(value);
			valueObject = new Integer(val);
			valueType = Integer.TYPE;
		    } catch (NumberFormatException ex) {
			throw new AppException("Invalid int: " + value);
		    }
		} else if (tag.equals("double")) {
		    String value = echild.getAttribute("value");
		    try {
			double val = Double.parseDouble(value);
			valueObject = new Double(val);
			valueType = Double.TYPE;
		    } catch (NumberFormatException ex) {
			throw new AppException("Invalid double: " + value);
		    }
		} else if (tag.equals("boolean")) {
		    String value = echild.getAttribute("value");
		    valueObject = new Boolean(value);
		    valueType = Boolean.TYPE;
		} else if (tag.equals("set")) {
		    valueObject = new HashSet();
		    configureBean(valueObject, echild, refs);
		    valueType = valueObject.getClass();
		} else if (tag.equals("list")) {
		    valueObject = new ArrayList();
		    configureBean(valueObject, echild, refs);
		    valueType = valueObject.getClass();
		} else if (tag.equals("map")) {
		    valueObject = new SmallMap();
		    configureBean(valueObject, echild, refs);
		    valueType = valueObject.getClass();
		} else if (tag.equals("annotation")) {
		    valueObject = new SmallAnnotation();
		    configureBean(valueObject, echild, refs);
		    valueType = valueObject.getClass();
		} else {
		    throw new AppException("Unknown element `" + tag + "' in XML-bean");
		}

		
		if (name != null && name.length() > 0) {
		    setProp(clazz, bean, name, valueObject, valueType);
		} else {
		    if (bean instanceof Collection) {
		        ((Collection) bean).add(valueObject);
		    } else {
		        throw new AppException("Anonymous beans are only allowed as children of Collections");
		    }
		}
		
	    }
	    child = child.getNextSibling();
	}
    }

    private void setProp(Class clazz, Object bean, String prop, Object value, Class ourType) 
        throws AppException
    {
		BeanInfo bi = null;
	
		try {
		    bi = Introspector.getBeanInfo(clazz);
		} catch (IntrospectionException ex) {
		    throw new AppException("Couldn't introspect class " + bean.getClass().getName());
		}
		PropertyDescriptor[] descs = bi.getPropertyDescriptors();
		for (int i = 0; i < descs.length; ++i) {
		    if (descs[i].getName().equals(prop)) {
			PropertyDescriptor desc = descs[i];
			if (! desc.getPropertyType().isAssignableFrom(ourType)) {
			    throw new AppException("Property " + prop + " is not assignable from " + ourType.getName());
			}
			Object[] obj = new Object[1];
			obj[0] = value;
			try {
			    desc.getWriteMethod().invoke(bean, obj);
			} catch (Exception ex) {
			    throw new AppException("Invocation failed, could not invoke " + prop + " " + value);
			}
			return;
		    }
		}
		if (bean instanceof Map) {
		    ((Map) bean).put(prop, value);
		} else if (bean instanceof Annotation) {
		    try {
		        ((Annotation) bean).setProperty(prop, value);
		    } catch (ChangeVetoException cve) {
		        throw new AppException("Unexpected veto updating Annotation");
		    }
		} else {
		    throw new AppException("Couldn't find property " + prop + " in class " + clazz.getName());
		}
    }
}
 
