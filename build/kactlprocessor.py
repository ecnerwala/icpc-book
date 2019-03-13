#!/usr/bin/env python2
# encoding: utf-8

# Source code preprocessor for KACTL building process.
# Written by Håkan Terelius, 2008-11-24

from __future__ import print_function
import sys
import getopt
import subprocess


def escape(input):
    input = input.replace('<', r'\ensuremath{<}')
    input = input.replace('>', r'\ensuremath{>}')
    return input

def pathescape(input):
    input = input.replace('\\', r'\\')
    input = input.replace('_', r'\_')
    input = escape(input)
    return input

def codeescape(input):
    input = input.replace('_', r'\_')
    input = input.replace('\n', '\\\\\n')
    input = input.replace('{', r'\{')
    input = input.replace('}', r'\}')
    input = input.replace('^', r'\ensuremath{\hat{\;}}')
    input = escape(input)
    return input

def ordoescape(input, esc=True):
    if esc:
        input = escape(input)
    start = input.find("O(")
    if start >= 0:
        bracketcount = 1
        end = start+1
        while end+1<len(input) and bracketcount>0:
            end = end + 1
            if input[end] == '(':
                bracketcount = bracketcount + 1
            elif input[end] == ')':
                bracketcount = bracketcount - 1
        if bracketcount == 0:
            return r"%s\bigo{%s}%s" % (input[:start], input[start+2:end], ordoescape(input[end+1:], False))
    return input

def addref(caption, outstream):
    caption = pathescape(caption).strip()
    print(r"\kactlref{%s}" % caption, file=outstream)
    with open('header.tmp', 'a') as f:
        f.write(caption + "\n")

def processwithcomments(caption, instream, outstream, listingslang = None):
    knowncommands = ['Author', 'Date', 'Description', 'Source', 'Time', 'Memory', 'License', 'Status', 'Usage']
    requiredcommands = ['Author', 'Description']
    includelist = []
    error = ""
    warning = ""
    # Read lines from source file
    try:
        lines = instream.readlines()
    except:
        error = "Could not read source."
    nlines = list()
    cur_hash = None
    for line in lines:
        had_comment = "///" in line
        # Remove /// comments
        if had_comment:
            line, tail = line.split("///",1)
            tail = tail.strip()
        line = line.rstrip()
        # Remove '#pragma once' and 'using namespace std;' lines
        if line == "#pragma once" or line == "using namespace std;":
            continue
        if line.endswith("/** exclude-line */"):
            continue
        if had_comment and not line and tail not in ["start-hash", "end-hash"]:
            continue
        # Check includes
        include = isinclude(line)
        if include is not None:
            includelist.append(include)
            continue

        if had_comment and tail == "start-hash":
            if line: line += ' '
            line += "// start-hash"
            cur_hash = []

        if cur_hash is not None: cur_hash.append(line)

        if had_comment and tail == "end-hash":
            script = 'hash' if listingslang else 'hash-cpp'
            cur_hash = '\n'.join(cur_hash)
            p = subprocess.Popen(['sh', '../content/contest/%s.sh' % script], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            hsh, _ = p.communicate(cur_hash)
            if line: line += ' '
            line += "// %s = %s" % (script, hsh.split(None,1)[0])
            cur_hash = None

        nlines.append(line)
    # Remove and process /** */ comments
    source = '\n'.join(nlines)
    nsource = ''
    start = source.find("/**")
    end = 0
    commands = {}
    while start >= 0 and not error:
        nsource = nsource.rstrip() + source[end:start]
        end = source.find("*/", start)
        if end<start:
            error = "Invalid /** */ comments."
            break
        comment = source[start+3:end].strip()
        end = end + 2
        start = source.find("/**",end)

        commentlines = comment.split('\n')
        command = None
        value = ""
        for cline in commentlines:
            allow_command = False
            cline = cline.strip()
            if cline.startswith('*'):
                cline = cline[1:].strip()
                allow_command = True
            ind = cline.find(':')
            if allow_command and ind != -1 and ' ' not in cline[:ind]:
                if command:
                    if command not in knowncommands:
                        error = error + "Unknown command: " + command + ". "
                    commands[command] = value.lstrip()
                command = cline[:ind]
                value = cline[ind+1:].strip()
            else:
                value = value + '\n' + cline
        if command:
            if command not in knowncommands:
                error = error + "Unknown command: " + command + ". "
            commands[command] = value.lstrip()
    for rcommand in sorted(set(requiredcommands) - set(commands)):
        error = error + "Missing command: " + rcommand + ". "
    if end>=0:
        nsource = nsource.rstrip() + source[end:]
    nsource = nsource.strip()

    # Produce output
    out = []
    if warning:
        out.append(r"\kactlwarning{%s: %s}" % (caption, warning))
    if error:
        out.append(r"\kactlerror{%s: %s}" % (caption, error))
    else:
        addref(caption, outstream)
        if commands.get("Description"):
            out.append(r"\defdescription{%s}" % escape(commands["Description"]))
        if commands.get("Usage"):
            out.append(r"\defusage{%s}" % codeescape(commands["Usage"]))
        if commands.get("Time"):
            out.append(r"\deftime{%s}" % ordoescape(commands["Time"]))
        if commands.get("Memory"):
            out.append(r"\defmemory{%s}" % ordoescape(commands["Memory"]))
        if includelist:
            out.append(r"\leftcaption{%s}" % pathescape(", ".join(includelist)))
        if nsource:
            out.append(r"\rightcaption{%d lines}" % len(nsource.split("\n")))
        langstr = ", language="+listingslang if listingslang else ""
        out.append(r"\begin{lstlisting}[caption={%s}%s]" % (pathescape(caption), langstr))
        out.append(nsource)
        out.append(r"\end{lstlisting}")

    for line in out:
        print(line, file=outstream)

def processraw(caption, instream, outstream, listingslang = 'raw'):
    try:
        source = instream.read().strip()
        addref(caption, outstream)
        print(r"\rightcaption{%d lines}" % len(source.split("\n")), file=outstream)
        print(r"\begin{lstlisting}[language=%s,caption={%s}]" % (listingslang, pathescape(caption)), file=outstream)
        print(source, file=outstream)
        print(r"\end{lstlisting}", file=outstream)
    except:
        print("\kactlerror{Could not read source.}", file=outstream)

def isinclude(line):
    line = line.strip()
    if line.startswith("#include") and not line.endswith("/** keep-include */"):
        return line[8:].strip()
    return None

def getlang(input):
    return input.rsplit('.',1)[-1]

def getfilename(input):
    return input.rsplit('/',1)[-1]

def print_header(data, outstream):
    parts = data.split('|')
    until = parts[0].strip() or parts[1].strip()
    if not until:
        # Nothing on this page, skip it.
        return
    with open('header.tmp') as f:
        lines = [x.strip() for x in f.readlines()]
    if until not in lines:
        # Nothing new on the page.
        return

    ind = lines.index(until) + 1
    def adjust(name):
        return name if name.startswith('.') else name.split('.')[0]
    output = r"\enspace{}".join(map(adjust, lines[:ind]))
    print(output, file=outstream)
    with open('header.tmp', 'w') as f:
        for line in lines[ind:]:
            f.write(line + "\n")

def main():
    language = None
    caption = None
    instream = sys.stdin
    outstream = sys.stdout
    print_header_value = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:l:c:", ["help", "output=", "input=", "language=", "caption=", "print-header="])
        for option, value in opts:
            if option in ("-h", "--help"):
                print("This is the help section for this program")
                print()
                print("Available commands are:")
                print("\t -o --output")
                print("\t -h --help")
                print("\t -i --input")
                print("\t -l --language")
                print("\t --print-header")
                return
            if option in ("-o", "--output"):
                outstream = open(value, "w")
            if option in ("-i", "--input"):
                instream = open(value)
                if language == None:
                    language = getlang(value)
                if caption == None:
                    caption = getfilename(value)
            if option in ("-l", "--language"):
                language = value
            if option in ("-c", "--caption"):
                caption = value
            if option == "--print-header":
                print_header_value = value
        if print_header_value is not None:
            print_header(print_header_value, outstream)
            return
        print(" * \x1b[1m{}\x1b[0m".format(caption))
        if language == "cpp" or language == "cc" or language == "c" or language == "h" or language == "hpp":
            processwithcomments(caption, instream, outstream)
        elif language == "java":
            processwithcomments(caption, instream, outstream, 'Java')
        elif language == "ps":
            processraw(caption, instream, outstream) # PostScript was added in listings v1.4
        elif language == "raw":
            processraw(caption, instream, outstream)
        elif language == "rawcpp":
            processraw(caption, instream, outstream, 'C++')
        elif language == "sh":
            processraw(caption, instream, outstream, 'bash')
        elif language == "py":
            processraw(caption, instream, outstream, 'Python')
        else:
            raise ValueError("Unkown language: " + str(language))
    except (ValueError, getopt.GetoptError, IOError) as err:
        print(str(err), file=sys.stderr)
        print("\t for help use --help", file=sys.stderr)
        return 2

if __name__ == "__main__":
    exit(main())
