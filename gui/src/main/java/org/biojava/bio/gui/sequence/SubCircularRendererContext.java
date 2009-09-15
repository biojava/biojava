package org.biojava.bio.gui.sequence;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.SymbolList;

/**
 * A renderer context that allows some or all properties of another context to
 * be over-ridden.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class SubCircularRendererContext
implements CircularRendererContext {
  private final CircularRendererContext delegate;
  private final SymbolList symbols;
  private final FeatureHolder features;
  private final double radius;

  /**
   * Create a new sub context.
   *
   * <p>
   * Supply the real values for symbols, features and radius if you want this
   * context to mask the values of the parent context. Otherwise, provide the
   * default values.
   * </p>
   *
   * @param delegate  the original context to wrap
   * @param symbols   the SymbolList to return for getSymbols(), or null
   * @param features  the FeatureHolder to return for getFeatures(), or null
   * @param radius    the radius to return for getRadius(), or NaN
   * @throws NullPointerException  if delegate is null
   */
  public SubCircularRendererContext(
          CircularRendererContext delegate,
          SymbolList symbols,
          FeatureHolder features,
          double radius)
  {
    if(delegate == null) {
      throw new NullPointerException("Delegate can not be null");
    }

    this.delegate = delegate;
    this.symbols = symbols;
    this.features = features;
    this.radius = radius;
  }

  public double getOffset() {
    return delegate.getOffset();
  }

  public double getAngle(int indx) {
    return delegate.getAngle(indx);
  }

  public int getIndex(double angle) {
    return 0;
  }

  public double getRadius() {
    if(Double.isNaN(radius)) {
      return delegate.getRadius();
    } else {
      return radius;
    }
  }

  public SymbolList getSymbols() {
    if(symbols == null) {
      return delegate.getSymbols();
    } else {
      return symbols;
    }
  }

  public FeatureHolder getFeatures() {
    if(features == null) {
      return delegate.getFeatures();
    } else {
      return features;
    }
  }
}
