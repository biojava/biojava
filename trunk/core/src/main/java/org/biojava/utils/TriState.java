

package org.biojava.utils;

/**
 * Class that represents the tristate values possible in
 * a logical operation: true, false and indeterminate.
 *
 * @author David Huen
 */
public final class TriState
{
    public static final TriState TRUE = new TriState("TRUE");
    public static final TriState FALSE = new TriState("FALSE");
    public static final TriState INDETERMINATE = new TriState("INDETERMINATE");


    private TriState(String description)
    {
    }
}

