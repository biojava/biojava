



package org.biojava.utils.regex;

import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

class PatternChecker
{
    private static class Tokenizer
    {
        private String packedTxt;
        private int ptr = 0;

        final static int EOL = -1;
        final static int SYMBOL_TOKEN = 0;
        final static int NUMERIC = 1;
        final static int LEFT_BRACE = 2;
        final static int RIGHT_BRACE = 3;
        final static int COMMA = 4;
        final static int LEFT_BRACKET = 5;
        final static int RIGHT_BRACKET = 6;
        final static int PLUS = 7;
        final static int ASTERISK = 8;
        final static int VERT_BAR = 9;
        final static int LEFT_SQBRACKET = 10;
        final static int RIGHT_SQBRACKET = 11;
        final static int Q_MARK = 12;
        final static int CARET = 13;
        final static int DOLLAR = 14;
        final static int DOT = 15;
        final static int UNKNOWN = 999;

        private Tokenizer(String target)
        {
            packedTxt = pack(target);
        }

        private char getToken()
            throws IndexOutOfBoundsException
        {
            if (hasNext())
                 return packedTxt.charAt(ptr++);
            else
                throw new IndexOutOfBoundsException("text length: " + packedTxt.length() + " index: " + ptr);
        }

        private int nextTokenType()
        {
            if (!hasNext()) return EOL;

            char nextCh = packedTxt.charAt(ptr);

            // symbol tokens assumed to be alphas.
            if (Character.isLetter(nextCh))
                return SYMBOL_TOKEN;
            if (Character.isDigit(nextCh))
                return NUMERIC;
            // now check for specific chars
            if (nextCh == '.')
                return DOT;
            if (nextCh == '{')
                return LEFT_BRACE;
            if (nextCh == '}')
                return RIGHT_BRACE;
            if (nextCh == ',')
                return COMMA;
            if (nextCh == '(')
                return LEFT_BRACKET;
            if (nextCh == ')')
                return RIGHT_BRACKET;
            if (nextCh == '+')
                return PLUS;
            if (nextCh == '*')
                return ASTERISK;
            if (nextCh == '|')
                return VERT_BAR;
            if (nextCh == '[')
                return LEFT_SQBRACKET;
            if (nextCh == ']')
                return RIGHT_SQBRACKET;
            if (nextCh == '?')
                return Q_MARK;
            if (nextCh == '^')
                return CARET;
            if (nextCh == '$')
                return DOLLAR;
            return UNKNOWN;
        }

        private boolean hasNext()
        {
            return ptr < packedTxt.length();
        }

        /**
         * produces a version of the String with whitespace removed.
         */
        private String pack(String source)
        {
            StringBuffer packedString = new StringBuffer();

            for (int i=0; i < source.length(); i++) {
                char currCh;
                if (!Character.isWhitespace(currCh = source.charAt(i))) {
                    packedString.append(currCh);
                }
            }

            return packedString.toString();
        }
    }

    private SymbolTokenization sToke = null;
    private Tokenizer toke = null;
    private StringBuffer output = null;
    private FiniteAlphabet alfa = null;

    PatternChecker(FiniteAlphabet alfa)
    {
        this.alfa = alfa;
    }

    String parse(String patternString)
        throws RegexException, IllegalAlphabetException
    {
        output = new StringBuffer();

        try {
            sToke = alfa.getTokenization("token");
            if (sToke.getTokenType() != SymbolTokenization.CHARACTER)
                throw new IllegalAlphabetException("This alphabet does not have a character tokenization.");
        }
        catch (BioException be) {
            throw new IllegalAlphabetException(be, "Attempt to get a tokenization for this alphabet failed.");
        }

        toke = new Tokenizer(patternString);
        parse();
        return output.toString();
    }

    void parse()
        throws RegexException
    {
        boolean gotContent = false;
        while (toke.hasNext()) {
            int tokenType = toke.nextTokenType();
            switch (tokenType) {
                case Tokenizer.LEFT_BRACKET:
                    //System.out.println("LEFT_BRACKET");
                    gotContent = true;
                    output.append(toke.getToken());    // consume left bracket
                    parse();
                    if (toke.nextTokenType() == Tokenizer.RIGHT_BRACKET) output.append(toke.getToken());
                    else throw new RegexException("syntax error in regex: right bracket expected.");
                    parseQuantifier();
                    break;
                case Tokenizer.SYMBOL_TOKEN:
                    //System.out.println("SYMBOL_TOKEN");
                    gotContent = true;
                    try {
                        parseSymbol();
                    }
                    catch (IllegalSymbolException ise) { throw new RegexException(ise); }
                    parseQuantifier();
                    break;
                case Tokenizer.DOT:
                    gotContent = true;
                    output.append(toke.getToken());
                    parseQuantifier();
                    break;
                case Tokenizer.LEFT_SQBRACKET: // pick up variant symbols
                    //System.out.println("LEFT_SQBRACKET");
                    gotContent = true;
                    parseVariantSymbols();
                    parseQuantifier();
                    break;
                case Tokenizer.VERT_BAR:
                    //System.out.println("VERT_BAR");
                    if (gotContent) {
                        gotContent = false;
                        output.append(toke.getToken());
                    }
                    else throw new RegexException("syntax error in specifying alternate pattern.");
                    break;
                case Tokenizer.EOL:
                    //System.out.println("EOL");
                    return;
                case Tokenizer.RIGHT_BRACKET:
                    //System.out.println("RIGHT_BRACKET");
                    if (gotContent)
                        return;
                    else throw new RegexException("unexpected right bracket.");
                default:
                    throw new RegexException("unexpected symbol " + (toke.getToken()));
            }

        }
    }

    void parseSymbol()
        throws IllegalSymbolException
    {
        // we need to parse the token to a symbol
        // check for ambiguity and make explicit
        // that ambiguity
        while (toke.nextTokenType() == Tokenizer.SYMBOL_TOKEN) {
            char token = toke.getToken();
            Symbol sym = sToke.parseToken(Character.toString(token));
    
            if (sym instanceof AtomicSymbol) {
                output.append(token);
            }
            else {
                // ambiguity
                // sort symbols according to name
                SortedSet sortedAmbigs = new TreeSet(
                    new Comparator () {
                        public int compare(Object o1, Object o2)
                        {
                            if ((o1 instanceof AtomicSymbol) &&
                                (o2 instanceof AtomicSymbol)) {
                                try {
                                    return (int) (sToke.tokenizeSymbol((AtomicSymbol) o1).charAt(0) 
                                        - sToke.tokenizeSymbol((AtomicSymbol) o2).charAt(0));
                                }
                                catch (IllegalSymbolException ise) {
                                    throw new AssertionError(ise);
                                }
                            }
                            else
                                throw new ClassCastException();
                        }
                    } 
                    );
                for (Iterator symI = ((FiniteAlphabet) sym.getMatches()).iterator();
                    symI.hasNext(); ) {
                    Symbol atomicSym = (Symbol) symI.next();
                    sortedAmbigs.add(atomicSym);
                }
    
                output.append("[");
                for (Iterator symI = sortedAmbigs.iterator();
                    symI.hasNext(); ) {
                    Symbol atomicSym = (Symbol) symI.next();
                    output.append(sToke.tokenizeSymbol(atomicSym));
                }
                output.append("]");
            }
        }
    }

    void parseQuantifier()
        throws RegexException
    {
        int tokenType = toke.nextTokenType();

        switch (tokenType) {
            case Tokenizer.LEFT_BRACE:
                getIterations();
                break;
            case Tokenizer.PLUS:
                output.append(toke.getToken());
                break;
            case Tokenizer.ASTERISK:
                output.append(toke.getToken());
                break;
            case Tokenizer.Q_MARK:
                output.append(toke.getToken());
                break;
            default:
        }

        // check for non-greediness
        if (toke.nextTokenType() == Tokenizer.Q_MARK) {
            output.append(toke.getToken());
        }
    }

    void getIterations()
        throws RegexException
    {
        // get left brace
        output.append(toke.getToken());

        // there can be one or two numerical arguments
        // {m,n} {m,} {,n} {m}
        int argCount = 0;
        boolean gotContent = false;

        while (toke.hasNext()) {
            int tokenType = toke.nextTokenType();

            switch (tokenType) {
                case Tokenizer.NUMERIC:
                    if (argCount > 1) throw new RegexException("too many arguments in quantifier");
                    // consume
                    while (toke.nextTokenType() == Tokenizer.NUMERIC) {
                        gotContent = true;
                        output.append(toke.getToken());
                    }
                    if (toke.nextTokenType() == Tokenizer.EOL)
                        throw new RegexException("syntax error: unexpected EOL");
                    break;
                case Tokenizer.COMMA:
                    if (argCount++ != 0) throw new RegexException("too many arguments in quantifier");
                    output.append(toke.getToken());
                    break;
                case Tokenizer.RIGHT_BRACE:
                    if (argCount++ > 1) throw new RegexException("too many arguments in quantifier");
                    if (!gotContent) throw new RegexException("no arguments were actually specified!");
                    output.append(toke.getToken());
                    return;
                default:
                    throw new RegexException("syntax error: unexpected symbol " + toke.getToken());
            }
        }
    }

    void parseVariantSymbols()
        throws RegexException
    {
        // consume left square bracket
        output.append(toke.getToken());

        while (toke.nextTokenType() == Tokenizer.SYMBOL_TOKEN) {

            char token = toke.getToken();

            Symbol sym;
            try {
                sym = sToke.parseToken(Character.toString(token));
            }
            catch (IllegalSymbolException ise) {
                throw new RegexException(ise);
            }

            if (sym instanceof AtomicSymbol) {
                output.append(token);
            }
            else {
                throw new RegexException("all variant symbols must be atomic.");
            }
        }

        // check and consume right bracket
        if (toke.nextTokenType() == Tokenizer.RIGHT_SQBRACKET) output.append(toke.getToken());
        else throw new RegexException("missing right square bracket while specifying variants. Encountered " + toke.getToken() + " instead.");
    }

}



 
