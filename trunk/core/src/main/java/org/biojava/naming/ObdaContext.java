package org.biojava.naming;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.NoSuchElementException;

import javax.naming.Binding;
import javax.naming.CompositeName;
import javax.naming.Context;
import javax.naming.InvalidNameException;
import javax.naming.Name;
import javax.naming.NameClassPair;
import javax.naming.NameNotFoundException;
import javax.naming.NameParser;
import javax.naming.NamingEnumeration;
import javax.naming.NamingException;
import javax.naming.NotContextException;
import javax.naming.OperationNotSupportedException;
import javax.naming.directory.Attributes;
import javax.naming.directory.BasicAttributes;
import javax.naming.directory.DirContext;
import javax.naming.directory.ModificationItem;
import javax.naming.directory.SearchControls;

/**
 *
 *
 * @author Matthew Pocock
 */
public class ObdaContext
        implements DirContext
{
  private ObdaContext parent;
  private String myAtomicName;
  private Hashtable bindings;
  private Hashtable env;
  private BasicAttributes attrs;

  protected Name getMyComponents(Name name) throws NamingException
  {
    if(name instanceof CompositeName) {
      if(name.size() > 1) {
        throw new InvalidNameException(
                name.toString() +
                " has more components than namespace can handle");
      }

      // Turn component that belongs to you into a compound name
      return ObdaUriParser.getInstance().parse(name.get(0));
    } else {
      // Already parsed
      return name;
    }
  }

  ObdaContext(ObdaContext parent,
              String myAtomicName,
              Hashtable bindings,
              Hashtable env,
              BasicAttributes attrs)
  {
    this.parent = parent;
    this.myAtomicName = myAtomicName;
    this.bindings = new Hashtable(bindings);
    this.env = new Hashtable(env);
    this.attrs = (BasicAttributes) attrs.clone();
  }

  Hashtable getBindings()
  {
    return bindings;
  }

  BasicAttributes getAttrs()
  {
    return attrs;
  }

  public String getNameInNamespace() throws NamingException
  {
    ObdaContext ancestor = parent;

    // No ancestor; at root of namespace
    if(ancestor == null) {
      return "";
    }

    Name name = ObdaUriParser.getInstance().parse("");
    name.add(myAtomicName);

    // Get the parent's names
    while(ancestor != null && ancestor.myAtomicName != null) {
      name.add(0, ancestor.myAtomicName);
      ancestor = ancestor.parent;
    }

    return name.toString();
  }

  public Name composeName(Name name, Name prefix) throws NamingException
  {
    Name result;

    // Both are compound names; compose using compound name rules
    if(!(name instanceof CompositeName) &&
            !(prefix instanceof CompositeName)) {
      result = (Name) (prefix.clone());
      result.addAll(name);
      return new CompositeName().add(result.toString());
    }

    // Simplistic implementation; do not support federation
    throw new OperationNotSupportedException(
            "Do not support composing composite names");
  }

  public Object lookup(Name name) throws NamingException
  {
    System.err.println("lookup: '" + name + "' for " + bindings);
    if(name.isEmpty()) {
      // Ask to look up this context itself; create and return
      // a new instance that has its own independent environment.
      System.err.println("Empty - return copy");
      return new ObdaContext(parent, myAtomicName, bindings, env, attrs);
    }

    // Extract the components that belong to this namespace
    Name nm = getMyComponents(name);
    System.err.println("My component is " + nm);

    // Find the object in the internal hash table
    String start = nm.get(0);
    Object answer = bindings.get(start);
    if(answer == null) {
      throw new NameNotFoundException(name + " not found");
    }

    if(nm.size() == 1) {
      return answer;
    } else {
      return ((Context) answer).lookup(nm.getSuffix(1));
    }
  }

  public NamingEnumeration list(Name name) throws NamingException
  {
    if(name.isEmpty()) {
      // Generate enumeration of context's contents
      return new ListOfNames(bindings.keys());
    }

    // Perhaps "name" names a context
    Object target = lookup(name);
    if(target instanceof Context) {
      try {
        return ((Context) target).list("");
      } finally {
        ((Context) target).close();
      }
    }
    throw new NotContextException(name + " cannot be listed");
  }

  public NamingEnumeration listBindings(Name name) throws NamingException
  {
    if(name.isEmpty()) {
      // Generate enumeration of context's contents
      return new ListOfBindings(bindings.keys());
    }

    // Perhaps "name" names a context
    Object target = lookup(name);
    if(target instanceof Context) {
      try {
        return ((Context) target).list("");
      } finally {
        ((Context) target).close();
      }
    }
    throw new NotContextException(name + " cannot be listed");
  }

  public void bind(Name name, Object obj) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void rebind(Name name, Object obj) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void unbind(Name name) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void rename(Name oldName, Name newName) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public Context createSubcontext(Name name) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void destroySubcontext(Name name) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public NameParser getNameParser(Name name) throws NamingException
  {
    // Do lookup to verify that the name exists
    Object obj = lookup(name);
    if(obj instanceof Context) {
      ((Context) obj).close();
    }
    return ObdaUriParser.getInstance();
  }

  public String composeName(String name, String prefix)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public Object addToEnvironment(String propName, Object propVal)
          throws NamingException
  {
    return env.put(propName, propVal);
  }

  public Object removeFromEnvironment(String propName)
          throws NamingException
  {
    return env.remove(propName);
  }

  public Hashtable getEnvironment() throws NamingException
  {
    return new Hashtable(env);
  }

  public void close() throws NamingException
  {
    // nothing to do
  }



  public NamingEnumeration list(String name) throws NamingException
  {
    return list(new CompositeName(name));
  }

  public NamingEnumeration listBindings(String name) throws NamingException
  {
    return listBindings(new CompositeName(name));
  }

  public Object lookup(String name) throws NamingException
  {
    return lookup(new CompositeName(name));
  }

  public Object lookupLink(Name name) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public Attributes getAttributes(Name name) throws NamingException
  {
    return getAttributes(name, null);  // Same as attrIds == null
  }

  public Attributes getAttributes(Name name, String[] attrIds)
          throws NamingException
  {
    if(name.isEmpty()) {
      // Ask for the attributes of this context
      return (Attributes) attrs.clone();
    }

    // Extract the components that belong to this namespace
    Name nm = getMyComponents(name);
    String atom = nm.get(0);
    Object inter = bindings.get(atom);

    if(nm.size() == 1) {
      // Atomic name; find object in the internal data structure
      if(inter == null) {
        throw new NameNotFoundException(name + " not found");
      }

      if(inter instanceof DirContext) {
        return ((DirContext) inter).getAttributes("", attrIds);
      } else {
        // Fetch the object's attributes from this context
        return (Attributes) attrs.clone();
      }
    } else {
      // Intermediate name; consume the name in this context
      // and then continue
      if(!(inter instanceof DirContext)) {
        throw new NotContextException(atom + " does not name a dircontext");
      }
      return ((DirContext) inter).getAttributes(nm.getSuffix(1), attrIds);
    }
  }

  public void modifyAttributes(Name name, int mod_op, Attributes attrs)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void modifyAttributes(Name name, ModificationItem[] mods)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void bind(Name name, Object obj, Attributes attrs)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void rebind(Name name, Object obj, Attributes attrs)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public DirContext createSubcontext(Name name, Attributes attrs)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public DirContext getSchema(Name name) throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public DirContext getSchemaClassDefinition(Name name)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public NamingEnumeration search(Name name, Attributes matchingAttrs)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public NamingEnumeration search(Name name,
                                  String filterExpr,
                                  Object[] filterArgs,
                                  SearchControls cons)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public NamingEnumeration search(Name name,
                                  Attributes matchingAttributes,
                                  String[] attributesToReturn)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public NamingEnumeration search(Name name,
                                  String filter,
                                  SearchControls cons)
          throws NamingException
  {
    throw new OperationNotSupportedException();
  }

  public void bind(String name, Object obj) throws NamingException
  {
    bind(new CompositeName(name), obj);
  }

  public void rebind(String name, Object obj) throws NamingException
  {
    rebind(new CompositeName(name), obj);
  }

  public void unbind(String name) throws NamingException
  {
    unbind(new CompositeName(name));
  }

  public void rename(String oldName, String newName) throws NamingException
  {
    rename(new CompositeName(oldName), new CompositeName(newName));
  }

  public Context createSubcontext(String name) throws NamingException
  {
    return createSubcontext(new CompositeName(name));
  }

  public void destroySubcontext(String name) throws NamingException
  {
    destroySubcontext(new CompositeName(name));
  }

  public NameParser getNameParser(String name) throws NamingException
  {
    return getNameParser(new CompositeName(name));
  }

  public Object lookupLink(String name) throws NamingException
  {
    return lookupLink(new CompositeName(name));
  }

  public Attributes getAttributes(String name) throws NamingException
  {
    return getAttributes(new CompositeName(name));
  }

  public Attributes getAttributes(String name, String[] attrIds)
          throws NamingException
  {
    return getAttributes(new CompositeName(name), attrIds);
  }

  public void modifyAttributes(String name, int mod_op, Attributes attrs)
          throws NamingException
  {
    modifyAttributes(new CompositeName(name), mod_op, attrs);
  }

  public void modifyAttributes(String name, ModificationItem[] mods)
          throws NamingException
  {
    modifyAttributes(new CompositeName(name), mods);
  }

  public void bind(String name, Object obj, Attributes attrs)
          throws NamingException
  {
    bind(new CompositeName(name), attrs);
  }

  public void rebind(String name, Object obj, Attributes attrs)
          throws NamingException
  {
    rebind(new CompositeName(name), obj, attrs);
  }

  public DirContext createSubcontext(String name, Attributes attrs)
          throws NamingException
  {
    return createSubcontext(new CompositeName(name), attrs);
  }

  public DirContext getSchema(String name) throws NamingException
  {
    return getSchema(new CompositeName(name));
  }

  public DirContext getSchemaClassDefinition(String name)
          throws NamingException
  {
    return getSchemaClassDefinition(new CompositeName(name));
  }

  public NamingEnumeration search(String name,
                                  Attributes matchingAttributes,
                                  String[] attributesToReturn)
          throws NamingException
  {
    return search(new CompositeName(name), matchingAttributes, attributesToReturn);
  }

  public NamingEnumeration search(String name,
                                  Attributes matchingAttributes)
          throws NamingException
  {
    return search(new CompositeName(name), matchingAttributes);
  }

  public NamingEnumeration search(String name,
                                  String filter,
                                  SearchControls cons)
          throws NamingException
  {
    return search(new CompositeName(name), filter, cons);
  }

  public NamingEnumeration search(String name,
                                  String filterExpr,
                                  Object[] filterArgs,
                                  SearchControls cons)
          throws NamingException
  {
    return search(new CompositeName(name), filterExpr, filterArgs, cons);
  }

  class ListOfNames implements NamingEnumeration {
    protected Enumeration names;

    ListOfNames(Enumeration names)
    {
      this.names = names;
    }

    public boolean hasMoreElements()
    {
      try {
        return hasMore();
      } catch (NamingException e) {
        return false;
      }
    }

    public boolean hasMore() throws NamingException
    {
      return names.hasMoreElements();
    }

    public Object next() throws NamingException
    {
      String name = (String) names.nextElement();
      String className = bindings.get(name).getClass().getName();
      return new NameClassPair(name, className);
    }

    public Object nextElement()
    {
      try {
        return next();
      } catch (NamingException e) {
        throw new NoSuchElementException(e.toString());
      }
    }

    public void close()
    {
    }
  }

  class ListOfBindings implements NamingEnumeration {
    protected Enumeration names;

    ListOfBindings(Enumeration names)
    {
      this.names = names;
    }

    public boolean hasMoreElements()
    {
      try {
        return hasMore();
      } catch (NamingException e) {
        return false;
      }
    }

    public boolean hasMore() throws NamingException
    {
      return names.hasMoreElements();
    }

    public Object next() throws NamingException
    {
      String name = (String) names.nextElement();
      Object bound = bindings.get(name);
      return new Binding(name, bound);
    }

    public Object nextElement()
    {
      try {
        return next();
      } catch (NamingException e) {
        throw new NoSuchElementException(e.toString());
      }
    }

    public void close()
    {
    }
  }
 }
