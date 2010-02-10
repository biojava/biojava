/**
 * The package contains Entity representations of BioJava classes.
 * The purpose of these entities is to allow simple serialization of BioJava data
 * using binary serialization for protocols that require this (eg RPC between 
 * Java application servers) as well as persistence mechanisms that require bean
 * like ojbects such as the Java Persistence Architechture (JPA) or the
 * Java API for XML Binding (JAXB). For this reason all objects in this package
 * should provide a parameterless public constructor and public get/set methods
 * for relevant fields.
 * <p>
 * Given the public nature of the constructors and the setters in these beans
 * these classes are not intended for direct use in general programming when 
 * using the BioJava v3 API. This is because it is possible to leave the bean in
 * and inconsitent state and they are <b>not thread safe</b> unless synchronization 
 * controlled externally (via synchornization blocks or via a application container). 
 * </p><p>
 * The Entities are intended to back other objects that a 
 * programer will interact with directly. For example <code>Foo.class</code> will be backed
 * by <code>FooEntity.class</code>. Generally interaction with Foo.class is to be prefered and
 * will often be more sensible as the entities typically provide no 'biological
 * behaivour'. Relevant behaivour should be provided by the wrapping class. It is best
 * to think of <code>Foo</code> as a view onto the data that is held in the 
 * <code>FooEntity</code>.  A good example is the sophisticated Symbol
 * behaivour that can represent biological logic about IUPAC ambiguity symbols.
 * For example a 'w' in a Biosequence represents an abiguity between 'a' and 't',
 * whereas a 'w' in BiosequenceEntity is simply a 'w' and nothing else.
 * </p><p>
 * The wrapper entity pattern is intended to allow for a lot of the advanced
 * behaivour in the original BioJava while also allowing use of modern transport
 * and persistence packages. This is achieved by peristing and transporting the 
 * entity without the wrapper and re-wrapping it at the other end.
 * </p><p>
 * Currently BioJava v3 uses annotated @Id fields to define
 * <code>equals(Object o)</code>. Consistent definition is critical to how
 * the object will behave when persisted to a database. In the case of:
 * <pre>
 * Foo f = ... initialize
 * Foo fo = ... initialize
 * boolean b = f.equals(fo);
 * </pre>
 * <code>b</code> would be true if both objects share the same value 
 * (or embeddable object) in the field that represents the primary key in the 
 * database <b>even</b> if all other fields are equal. This is desirable because
 * two entities representing the same DB record may be retreived from two different
 * sessions. Additionally these are the identity fields, so logically, they should map to
 * the concept of identity. Finally, searching a collection is made very simple
 * without requireing an iterator:
 * <pre>
 * Integer id = //code to initialize 
 * collection.contains(new Foo(id));
 * </pre>
 * By default BioJava v3 entities use <b>only</b> the primary key field for equality
 * If either record has <code>null</code> as the primary key value it is never equal
 * to another. When implementing <code>equals(Object o)</code> it is not advisable to perform
 * the test this.getClass() == o.getClass() because of the possibility of proxy
 * classes used in JPA. This can, however, lead to an issue with the 
 * <code>hashcode()</code> method.  Consider the following code:
 * <pre>
 * Foo foo = new Foo() //no primary key
 * HashSet set = new HashSet();
 * set.add(foo);
 * // code here to persist Foo and consequently generate it's PK
 * boolean b = set.contains(foo);
 * </pre>
 * Because only the PK is used for equality, then the PK is used in the hashcode.
 * This means that <code>b</code> is probably going to be false because
 * it would have been stored in a hash bucket using the old hashcode that will
 * now be different even though the set actually does contain a pointer to foo.
 * Although a potential deficiency it is unlikely to be a major problem for 
 * BioJava v3 developers because using entity backed objects is prefered to direct
 * interaction with entities. If you need to use entities directly then use hashed
 * collections with caution.
 * 
 * <p>Wrapper classes can either delegate it's equals call to the underlying 
 * entity or it can do something that is more biologically sensible 
 * (as PK values are typically not exposed in the wrapper). It is probably more
 * sensible for a wrapper to define it's own <code>equals</code> (and <code>haschode</code>
 * implementations due to the limitations of the default @Id based system 
 * described above. Especially the potential hashcode problems.
 * 
 * For example <code>FooSequence.class</code> might want to base
 * equality on the exact match of the DNA sequence it holds even though
 * <code>FooSequenceEntity.class</code> may only use the PK field. If delegation
 * is used (or not) it should be clearly documented.
 * <p>
 * 
 * </p>
 * @author Mark Schreiber
 */
package org.biojava3.db.biosql;

