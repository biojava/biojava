


package org.biojava.utils.regex;

import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;

/**
 * A class analogous to java.util.regex.Pattern but for SymbolLists.
 * @author David Huen
 * @since 1.4
 */
public class Pattern
{
    private FiniteAlphabet alfa;
    private java.util.regex.Pattern pattern;
    private String label;
    private String patternTxt;

    Pattern(String patternTxt, PatternChecker checker, FiniteAlphabet alfa)
        throws IllegalAlphabetException, RegexException
    {
        this.pattern = java.util.regex.Pattern.compile(checker.parse(patternTxt));
        this.patternTxt = patternTxt;
        this.alfa = alfa;
        label = patternTxt;
    }

    Pattern(String patternTxt, PatternChecker checker, FiniteAlphabet alfa, String label)
        throws IllegalAlphabetException, RegexException
    {
        this.pattern = java.util.regex.Pattern.compile(checker.parse(patternTxt));
        this.patternTxt = patternTxt;
        this.alfa = alfa;
        this.label = label;
    }

    /**
     * return the String label associated with this pattern.
     */
    public String getName() { return label; }

    /**
     * Creates a matcher that will match the given input against this pattern.
     * @param sl SymbolList against which match is to be made.
     * @return A new matcher for this pattern.
     */
    public org.biojava.utils.regex.Matcher matcher(SymbolList sl)
    {
        return new org.biojava.utils.regex.Matcher(this, sl);
    }

    /**
     * returns the Pattern to be matched as a String.
     */
    public String patternAsString()
    {
        return patternTxt;
    }

    /**
     * returns the actual String used to construct the regex with all
     * ambiguities expanded.
    //FIXME: do something about unicode strings and conversion back to something sensible.
     */
    public String patternExpanded()
    {
        return pattern.pattern();
    }

    /**
     * returns the java.util.regex.Pattern object that underlies this instance.
     */
    java.util.regex.Pattern getPattern() { return pattern; }

    public FiniteAlphabet getAlphabet() { return alfa; }
}

