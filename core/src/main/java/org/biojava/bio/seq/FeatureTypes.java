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
package org.biojava.bio.seq;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.SmallMap;

/**
 * Registry of known types of features.
 *
 * <p>Historically, the type of a feature within a sequence has been represented
 * by the class or interface of the feature, or by the type property. Different
 * databases use different names for types, and the same name can mean different
 * things depending on where it came from.
 * </p>
 *
 * <p>This class attempts to provide a framework for registering and managing
 * feature types. It includes the concept of a type having a FeatureFilter
 * schema that would match all examples of that type, and a URI for that type.
 * Groups of related types, such as all those defined in embl, can be packaged
 * up into a single namespace. One type can refer to any number of others,
 * indicating that they are direct parents. This lets us say things like an
 * Ensembl exon is a specific type of the general exon type, and also an
 * instance of the Ensembl feature type.</p>
 *
 * <p>All feature types are presumed to have URIs that are of the form:
 * <code>uri:biojava.org:types:${repository}/${type}</code> where
 * <code>${repository}</code> is something like embl or ensembl. Feature types
 * defined within biojava will be in the repository called "core".
 * <code>${type}</code> is the local type of the feature. This will be something
 * like exon, or repeat. The total URI must be unique. Different repositories
 * can (and are encouraged to) define types with the same names.</p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class FeatureTypes {
  private static final Map repositories;
  
  /** The standard prefix for all type URIs **/
  public static final String URI_PREFIX = "uri:biojava.org:types";
  
  static {
    repositories = new SmallMap();
  }
  
  /**
   * <p>Fetch a repository by name.</p>
   *
   * <p>Use this to find out what types a repository can provide.</p>
   *
   * @param name  the name of the Repository to fetch
   * @return the Repository with that name
   * @throws NoSuchElementException  if there is no repository by that name
   */
  public static Repository getRepository(String name)
  throws NoSuchElementException {
    Repository rep = (Repository) repositories;
    
    if(rep == null) {
      throw new NoSuchElementException("Could not find repository: " + name);
    }
    
    return rep;
  }
  
  /**
   * Find the names of all known repositories.
   *
   * @return a Set of all repository names
   */
  public static Set getRepositoryNames() {
    return repositories.keySet();
  }
  
  /**
   * Add a repository to FeatureTypes.
   *
   * @param repos  the Repository to add
   */
  public static void addRepository(Repository repos) {
    repositories.put(repos.getName(), repos);
  }
  
  /**
   * Remove a repository from FeaturTypes.
   *
   * @param repos  the Repository to remove
   */
  public static void removeRepository(Repository repos) {
    repositories.remove(repos.getName());
  }
  
  /**
   * <p>Get a Type by URI.</p>
   *
   * <p>This will attemt to resolve the URI to a type by first matching the
   * standard prefix, then finding a repository of the right name and finally
   * searching that with the type-name portion of the URI.</p>
   *
   * @param uri  the URI to resolve
   * @return a Type with that URI
   * @throws NoSuchElementException if the type could not be resolved
   */
  public static Type getType(String uri) {
    if(!uri.startsWith(URI_PREFIX)) {
      throw new NoSuchElementException(
        "All types start with: " + URI_PREFIX +
        " while processing " + uri
      );
    }
    
    String names = uri.substring(URI_PREFIX.length() + 1);
    int slash = uri.indexOf("/");
    String repName = names.substring(0, slash);
    String typeName = names.substring(slash + 1);
    
    Repository rep = getRepository(repName);
    return rep.getType(typeName);
  }
  
  /**
   * <p>Work out if one type is a sub-type of another.</p>
   *
   * <p>This is the transiative closure of Type.getParent(). It will return true
   * if superType is reachable by following getParent() calls, and false
   * otherwise.</p>
   *
   * @param subType  the Type of the potential sub type
   * @param superType  the Type of the potential super type
   * @return true if they are sub-super types, false otherwise
   */
  public static boolean isSubTypeOf(Type subType, Type superType) {
    Set parents = subType.getParents();
    
    if(parents.contains(superType.getURI())) {
      return true;
    }
    
    for(Iterator i = parents.iterator(); i.hasNext(); ) {
      String puri = (String) i.next();
      return isSubTypeOf(getType(puri), superType);
    }
    
    return false;
  }
  
  /**
   * A named collection of Types.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static interface Repository
  extends Annotatable {
    /**
     * <p>The name of this repository.</p>
     *
     * <p>This will be the ${repository} component of any URIs of types defined
     * here.</p>
     *
     * @return the name of the repository
     */
    String getName();
    
    /**
     * Get a set of all type names defined in this repository.
     *
     * @return a Set of Type names as Strings
     */
    Set getTypes();
    
    /**
     * Find the type for a name.
     *
     * @param name  the name of the Type
     * @return  the Type of that name
     * @throws NoSuchElementException  if that type can not be found
     */
    Type getType(String name)
    throws NoSuchElementException;
  }
  
  /**
   * A type of feature.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static interface Type
  extends Annotatable {
    /**
     * <p>Get the schema for this type.</p>
     *
     * <p>The schema is represented as a FeatureFilter. This will almost
     * certainly be a complext filter using ands and ors to combine multiple
     * constraints. A particular type may chose to restrict any one of the
     * feature's properties, their allowed children and their allowed parents
     * in a feature hierachy, the type of the annotation associated with it and
     * anything else that can be expressed using a feature fitler.</p>
     *
     * <p>For a feature to actualy conform to this type, it must be acceptable
     * by the schema filter.</p>
     *
     * @return the schema FeatureFilter
     */
    FeatureFilter getSchema();
    
    /**
     * Get the name of this type.
     *
     * @return the Type name
     */
    String getName();
    
    /**
     * Get a set of URIs for parent types.
     *
     * @return a Set of all parent URIs
     */
    Set getParents();
    
    /**
     * <p>Get the URI for this type.</p>
     *
     * <p>The URI will be composed according to the rules defined in
     * FeatureTypes, being of the form
     * <code>uri:biojava.org:types:${repository}/${type}</code>.</p>
     *
     * @return the URI for this type
     */
    String getURI();
  }
  
  /**
   * A simple implementation of a Repository.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static class RepositoryImpl
  extends AbstractChangeable
  implements Repository {
    private final String name;
    private final Map types;
    private final Annotation ann;
    
    /**
     * Create a named repository.
     *
     * @param name the name of this repository
     */
    public RepositoryImpl(String name) {
      this.name = name;
      types = new HashMap();
      this.ann = new SimpleAnnotation();
    }
    
    protected ChangeSupport getChangeSupport(ChangeType ct) {
      ChangeSupport cs = super.getChangeSupport(ct);
      ann.addChangeListener(new Annotatable.AnnotationForwarder(this, cs), Annotation.PROPERTY);
      return cs;
    }
    
    public Annotation getAnnotation() {
      return ann;
    }
    
    public String getName() {
      return name;
    }
    
    public Set getTypes() {
      return types.keySet();
    }
    
    public Type getType(String name)
    throws NoSuchElementException {
      Type type = (Type) types.get(name);
      if(type == null) {
        throw new NoSuchElementException(
          "Could not find type " + name +
          " in repository " + getName()
        );
      }
      return type;
    }
    
    /**
     * Create a new type in this repository.
     *
     * @param name  the Type name
     * @param schema  the FeatureFilter defining the type
     * @param parents  the Set (possibly empty) of parent URIs
     */
    public Type createType(
      String name,
      FeatureFilter schema,
      Set parents
    ) {
      Type type = new TypeImpl(name, schema, parents);
      types.put(name, type);
      return type;
    }
    
    private class TypeImpl
    extends AbstractChangeable
    implements Type {
      private Annotation ann;
      private String name;
      private FeatureFilter schema;
      private Set parents;
      
      public TypeImpl(
        String name,
        FeatureFilter schema,
        Set parents
      ) {
        this.name = name;
        this.schema = schema;
        this.parents = parents;
        this.ann = new SimpleAnnotation();
      }
      
      protected ChangeSupport getChangeSupport(ChangeType ct) {
        ChangeSupport cs = super.getChangeSupport(ct);
        ann.addChangeListener(new Annotatable.AnnotationForwarder(this, cs), Annotation.PROPERTY);
        return cs;
      }
      
      public Annotation getAnnotation() {
        return ann;
      }
      
      public String getName() {
        return name;
      }
      
      public FeatureFilter getSchema() {
        return schema;
      }
      
      public Set getParents() {
        return parents;
      }
      
      public String getURI() {
        return URI_PREFIX + "/" + RepositoryImpl.this.getName() + "/" + name;
      }
    }
  }
}
