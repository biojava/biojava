// $Id: CommandLineArguments.java,v 1.16 2009/10/26 23:29:40 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public final class CommandLineArguments {

    private final static String OPTIONS_PREFIX          = "-";
    private final static String EXTENDED_OPTIONS_PREFIX = "--";
    private final static String OPTIONS_SEPARATOR       = "=";
    private Map<String, String> _options;
    private Map<String, String> _extended_options;
    private List<String>        _names;
    private String              _command_line_str;

    public CommandLineArguments( final String[] args ) throws IOException {
        init();
        parseCommandLineArguments( args );
    }

    private Map<String, String> getAllOptions() {
        final Map<String, String> o = new HashMap<String, String>();
        o.putAll( getOptionsList() );
        o.putAll( getExtendedOptionsList() );
        return o;
    }

    public String getCommandLineArgsAsString() {
        return _command_line_str;
    }

    private Map<String, String> getExtendedOptionsList() {
        return _extended_options;
    }

    public File getFile( final int i ) {
        return new File( getNames()[ i ] );
    }

    public String getName( final int i ) {
        return getNames()[ i ];
    }

    public String[] getNames() {
        final String[] a = new String[ getNamesList().size() ];
        return getNamesList().toArray( a );
    }

    private List<String> getNamesList() {
        return _names;
    }

    public int getNumberOfNames() {
        return getNames().length;
    }

    private Map<String, String> getOptionsList() {
        return _options;
    }

    public String getOptionValue( final String option_name ) throws IllegalArgumentException {
        final Map<String, String> o = getAllOptions();
        if ( o.containsKey( option_name ) ) {
            final String value = o.get( option_name );
            if ( !ForesterUtil.isEmpty( value ) ) {
                return value;
            }
            else {
                throw new IllegalArgumentException( "value for \"" + option_name + "\" is not set" );
            }
        }
        else {
            throw new IllegalArgumentException( "option \"" + option_name + "\" is not set" );
        }
    }

    /**
     * Removes quotes
     * 
     */
    public String getOptionValueAsCleanString( final String option_name ) throws IllegalArgumentException {
        return getOptionValue( option_name ).replaceAll( "\"", "" ).replaceAll( "\'", "" );
    }

    public double getOptionValueAsDouble( final String option_name ) throws IOException {
        double d = -Double.MAX_VALUE;
        try {
            d = new Double( getOptionValue( option_name ) ).doubleValue();
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "value for option \"" + option_name + "\" is expected to be of type double" );
        }
        return d;
    }

    public int getOptionValueAsInt( final String option_name ) throws IOException {
        int i = Integer.MIN_VALUE;
        try {
            i = new Integer( getOptionValue( option_name ) ).intValue();
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "value for option \"" + option_name + "\" is expected to be of type integer" );
        }
        return i;
    }

    public long getOptionValueAsLong( final String option_name ) throws IOException {
        long l = Long.MIN_VALUE;
        try {
            l = new Long( getOptionValue( option_name ) ).longValue();
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "value for option \"" + option_name + "\" is expected to be of type long" );
        }
        return l;
    }

    private void init() {
        _options = new HashMap<String, String>();
        _extended_options = new HashMap<String, String>();
        _names = new ArrayList<String>();
        _command_line_str = "";
    }

    public boolean isOptionHasAValue( final String option_name ) {
        final Map<String, String> o = getAllOptions();
        if ( o.containsKey( option_name ) ) {
            final String value = o.get( option_name );
            return ( !ForesterUtil.isEmpty( value ) );
        }
        else {
            throw new IllegalArgumentException( "option \"" + option_name + "\" is not set" );
        }
    }

    public boolean isOptionSet( final String option_name ) {
        final Map<String, String> o = getAllOptions();
        return ( o.containsKey( option_name ) );
    }

    public boolean isOptionValueSet( final String option_name ) throws IllegalArgumentException {
        final Map<String, String> o = getAllOptions();
        if ( o.containsKey( option_name ) ) {
            return !( ForesterUtil.isEmpty( o.get( option_name ) ) );
        }
        else {
            throw new IllegalArgumentException( "option \"" + option_name + "\" is not set" );
        }
    }

    private void parseCommandLineArguments( final String[] args ) throws IOException {
        for( int i = 0; i < args.length; ++i ) {
            final String arg = args[ i ].trim();
            _command_line_str += arg;
            if ( i < args.length - 1 ) {
                _command_line_str += " ";
            }
            if ( arg.startsWith( CommandLineArguments.EXTENDED_OPTIONS_PREFIX ) ) {
                parseOption( arg.substring( CommandLineArguments.EXTENDED_OPTIONS_PREFIX.length() ),
                             getExtendedOptionsList() );
            }
            else if ( arg.startsWith( CommandLineArguments.OPTIONS_PREFIX ) ) {
                parseOption( arg.substring( CommandLineArguments.OPTIONS_PREFIX.length() ), getOptionsList() );
            }
            else {
                getNamesList().add( arg );
            }
        }
    }

    private void parseOption( final String option, final Map<String, String> options_map ) throws IOException {
        final int sep_index = option.indexOf( CommandLineArguments.OPTIONS_SEPARATOR );
        if ( sep_index < 1 ) {
            if ( ForesterUtil.isEmpty( option ) ) {
                throw new IOException( "attempt to set option with an empty name" );
            }
            if ( getAllOptions().containsKey( option ) ) {
                throw new IOException( "attempt to set option \"" + option + "\" mutiple times" );
            }
            options_map.put( option, null );
        }
        else {
            final String key = option.substring( 0, sep_index );
            final String value = option.substring( sep_index + 1 );
            if ( ForesterUtil.isEmpty( key ) ) {
                throw new IllegalArgumentException( "attempt to set option with an empty name" );
            }
            //  if ( ForesterUtil.isEmpty( value ) ) {
            //      throw new IllegalArgumentException( "attempt to set option with an empty value" );
            //  }
            if ( getAllOptions().containsKey( key ) ) {
                throw new IllegalArgumentException( "attempt to set option \"" + key + "\" mutiple times [" + option
                        + "]" );
            }
            options_map.put( key, value );
        }
    }

    public List<String> validateAllowedOptions( final List<String> allowed_options ) {
        final Map<String, String> options = getAllOptions();
        final List<String> dissallowed = new ArrayList<String>();
        for( final String o : options.keySet() ) {
            if ( !allowed_options.contains( o ) ) {
                dissallowed.add( o );
            }
        }
        return dissallowed;
    }

    public String validateAllowedOptionsAsString( final List<String> allowed_options ) {
        final List<String> dissallowed = validateAllowedOptions( allowed_options );
        String dissallowed_string = "";
        for( final Iterator<String> iter = dissallowed.iterator(); iter.hasNext(); ) {
            dissallowed_string += "\"" + iter.next();
            if ( iter.hasNext() ) {
                dissallowed_string += "\", ";
            }
            else {
                dissallowed_string += "\"";
            }
        }
        return dissallowed_string;
    }

    public List<String> validateMandatoryOptions( final List<String> mandatory_options ) {
        final Map<String, String> options = getAllOptions();
        final List<String> missing = new ArrayList<String>();
        for( final String string : mandatory_options ) {
            final String ma = string;
            if ( !options.containsKey( ma ) ) {
                missing.add( ma );
            }
        }
        return missing;
    }

    public String validateMandatoryOptionsAsString( final List<String> mandatory_options ) {
        final List<String> missing = validateMandatoryOptions( mandatory_options );
        String missing_string = "";
        for( final Iterator<String> iter = missing.iterator(); iter.hasNext(); ) {
            missing_string += "\"" + iter.next();
            if ( iter.hasNext() ) {
                missing_string += "\", ";
            }
            else {
                missing_string += "\"";
            }
        }
        return missing_string;
    }
}
