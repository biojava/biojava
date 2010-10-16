package org.biojava.utils.walker;

/**
 * Things that will be shown filters.
 *
 * <p>
 * Visitors will be shown filters by Walker instances. The walker will take care
 * of traversing filters like And, Or and ByParent that wrap other filters.
 * </p>
 *
 * <p>
 * Visitor implementations must be public or package-private. This is because
 * the walker implementation needs to bind to methods directly, and not via a
 * publicly accessible interface.
 * </p>
 * 
 * <p>
 * The simplest form of a
 * visitor has the single method void featureFilter(FeatureFilter ff). In this
 * case, the visitor will have this one method called for all filters.
 * </p>
 *
 * <h2>Writing Handler Methods.</h2>
 *
 * <p>
 * Visitor Handler methods should all follow the same pattern. This will allow
 * them to be recognised by introspection. For a feature filter class called
 * <code>Foo</code> or called <code>FeatureFilter.Foo</code>, the visitor
 * method should have the signature:
 * <pre>public void foo(Foo)</pre>
 * or
 * <pre>public Bar foo(Foo)</pre>
 * It is an error to provide both handlers, but then your compiler will tell
 * you this in advance.
 * </p>
 *
 * <p>
 * A given visitor should either provide void handlers, or always provide
 * handlers with the same return type. Do not mix these two. If you do, then
 * either the walker will not be created by WalkerFactory, or you will
 * experience strange run-time problems.
 * </p>
 *
 * <h2>How Handler Methods Will Be Chosen</h2>
 *
 * <p>
 * There may potentialy be several methods that could be used to handle a
 * feature filter. Using the above example, <code>Foo</code> obviously extends
 * <code>FeatureFilter</code>, so if there is a handler method for both classes,
 * there is an ambiguity. The tree walker will always call the handler method
 * that implements the most derived class. Only one handler will ever be
 * called for one filter.
 * </p>
 *
 * <p>
 * If there is no handler for a particular filter class, then no action will
 * be taken. In particular, if you supply a handler for just <code>Foo</code>,
 * and for no other filter types, even <code>FeatureFilter</code>, then only
 * <code>Foo</code> filters will be shown to the visitor. In most cases, it will
 * be wise to provide a handler for <code>FeatureFilter</code> even if you think
 * you handle all filter types, because this lets you throw a warning message.
 * </p>
 *
 * <h2>Handler Methods for <code>And</code> and Other Things</h2>
 *
 * <p>
 * Some feature filters, such as <code>FeatureFilter.And</code> wrap other
 * filters. In some circumstances, it is usefull to know about these, and
 * possibly to know about the results of visiting their children. The handler
 * for filters with child filters will always be invoked after all children
 * have been processed.
 * </p>
 *
 * <p>
 * In the case where all handlers return void, the handler for filters with
 * children will follow the normal pattern. So, the handler for
 * <code>FeatureFilter.And</code> will be
 * <pre>public void and(FeatureFilter.And)</pre>
 * In the case where all handlers return instances of some class, the handlers
 * for filters with children will contain one extra argument per child. So, the
 * handler in this case would be
 * <pre>public Bar and(FeatureFilter.And, Bar, Bar)</pre>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface Visitor {
}
