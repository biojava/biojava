# Development hints

Please use the provided formatters and templates to ensure a consistent code style.

Set up your development environment (assuming the use of Eclipse) as follows:

* import `biojava-formatter.xml` in `Preferences > Java > Code Style > Formatter`
* import `biojava-clenaup.xml` in `Preferences > Java > Code Style > Clean Up`
* import `biojava-templates.xml` in `Preferences > Java > Code Style > Code Templates`
* set `General > Workspace > Text` file encoding to `UTF-8`
* activate `Java > Editor > save Actions`: at least Format all lines, organize imports, the rest as you see fit, we suggest:
	* Convert control statement bodies to block
	* Add final modifier to private fields
	* Add final modifier to method parameters
	* Add final modifier to local variables
	* Remove unused imports
	* Add missing '@Override' annotations
	* Add missing '@Override' annotations to implementations of interface methods
	* Add missing '@Deprecated' annotations
	* Remove unnecessary casts
	* Remove unnecessary '$NON-NLS$' tags
	* Remove trailing white spaces on all lines
* configure XML formatting in `Preferences > XML > XML Files > Editor` 
	* Line width: `120`
	* Format comments (and join lines)
	* Insert whitespace before closing empty end-tags
	* Indent using `2` spaces
* install the Checkstyle Plugin and import `biojava-checkstyle.xml`
* (install the Findbugs Plugin)

Other useful eclipse plugins (find these via the Eclipse Marketplace):

* AnyEdit tools
* ObjectAid Class Diagram plugin
* Eclipse Mylyn Github Feature
