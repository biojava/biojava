package org.biojava.bio.symbol;

import junit.framework.TestCase;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Regression test for event forwarding bug found by James Worthington.
 *
 * @author Matthew Pocock
 */
public abstract class AbstractSymbolEventTest
        extends TestCase
{

  Annotation ann;
  Symbol sym;

  protected void setUp()
          throws Exception
  {
    Annotation initialAnno = new SmallAnnotation();
    sym = createSymbol(initialAnno);
    ann = sym.getAnnotation();
  }

  protected abstract Symbol createSymbol(Annotation ann) throws Exception;

  protected void doChange()
  {
    try {
      ann.setProperty("pigs", "dogs");
    } catch (ChangeVetoException e) {
      throw new AssertionFailure(e);
    }
  }

  public void testUnknown()
  {
    EventCounter uc = new EventCounter("unknown");
    sym.addChangeListener(uc, ChangeType.UNKNOWN);

    doChange();

    assertEquals("Prechange fired", 1, uc.getPre());
    assertEquals("Postchange fired", 1, uc.getPost());
  }

  public void testAnnotation()
  {
    EventCounter an = new EventCounter("annotation");
    sym.addChangeListener(an, Annotatable.ANNOTATION);
    doChange();

    assertEquals("Prechange fired", 1, an.getPre());
    assertEquals("Postchange fired", 1, an.getPost());
  }

  public void testProperty()
  {
    EventCounter pr = new EventCounter("property");
    sym.addChangeListener(pr, Annotation.PROPERTY);
    doChange();

    assertEquals("Prechange fired", 0, pr.getPre());
    assertEquals("Postchange fired", 0, pr.getPost());
  }

  public void testUnknownAnnotation()
  {
    EventCounter uc = new EventCounter("unknown");
    EventCounter an = new EventCounter("annotation");
    sym.addChangeListener(uc, ChangeType.UNKNOWN);
    sym.addChangeListener(an, Annotatable.ANNOTATION);

    doChange();

    assertEquals("Prechange fired", 1, uc.getPre());
    assertEquals("Postchange fired", 1, uc.getPost());
    assertEquals("Prechange fired", 1, an.getPre());
    assertEquals("Postchange fired", 1, an.getPost());
  }

  public void testUnknownProptery()
  {
    EventCounter uc = new EventCounter("unknown");
    EventCounter pr = new EventCounter("property");
    sym.addChangeListener(uc, ChangeType.UNKNOWN);
    sym.addChangeListener(pr, Annotation.PROPERTY);

    doChange();

    assertEquals("Prechange fired", 1, uc.getPre());
    assertEquals("Postchange fired", 1, uc.getPost());
    assertEquals("Prechange fired", 0, pr.getPre());
    assertEquals("Postchange fired", 0, pr.getPost());
  }

  public void testUnknownAnnotationProperty()
  {
    EventCounter uc = new EventCounter("unknown");
    EventCounter an = new EventCounter("annotation");
    EventCounter pr = new EventCounter("property");
    sym.addChangeListener(uc, ChangeType.UNKNOWN);
    sym.addChangeListener(an, Annotatable.ANNOTATION);
    sym.addChangeListener(pr, Annotation.PROPERTY);

    doChange();

    assertEquals("Prechange fired", 1, uc.getPre());
    assertEquals("Postchange fired", 1, uc.getPost());
    assertEquals("Prechange fired", 1, an.getPre());
    assertEquals("Postchange fired", 1, an.getPost());
    assertEquals("Prechange fired", 0, pr.getPre());
    assertEquals("Postchange fired", 0, pr.getPost());
  }

  public void testAnnotationUnknown()
  {
    EventCounter an = new EventCounter("annotation");
    EventCounter uc = new EventCounter("unknown");
    sym.addChangeListener(an, Annotatable.ANNOTATION);
    sym.addChangeListener(uc, ChangeType.UNKNOWN);

    doChange();

    assertEquals("Prechange fired", 1, uc.getPre());
    assertEquals("Postchange fired", 1, uc.getPost());
    assertEquals("Prechange fired", 1, an.getPre());
    assertEquals("Postchange fired", 1, an.getPost());
  }

  private class EventCounter
          implements ChangeListener
  {
    private final String name;

    private int pre = 0;
    private int post = 0;

    public EventCounter(String name)
    {
      this.name = name;
    }

    public void preChange(ChangeEvent cev)
            throws ChangeVetoException
    {
      pre++;
    }

    public void postChange(ChangeEvent cev)
    {
      post++;
    }

    public void clear()
    {
      pre = 0;
      post = 0;
    }

    public int getPre()
    {
      return pre;
    }

    public int getPost()
    {
      return post;
    }

    public String getName()
    {
      return name;
    }

    public String toString()
    {
      return name;
    }
  }
}