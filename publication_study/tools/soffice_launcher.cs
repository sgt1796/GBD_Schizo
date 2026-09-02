using System;
using System.Diagnostics;
using System.IO;
using System.Linq;

internal static class SofficeLauncher
{
    private static string Quote(string value)
    {
        return "\"" + value.Replace("\"", "\\\"") + "\"";
    }

    public static int Main(string[] args)
    {
        string programDirectory = Environment.GetEnvironmentVariable("CODEX_LIBREOFFICE_PROGRAM");
        if (String.IsNullOrWhiteSpace(programDirectory))
        {
            Console.Error.WriteLine("CODEX_LIBREOFFICE_PROGRAM must point to a LibreOffice program directory.");
            return 2;
        }
        string executable = Path.Combine(programDirectory, "soffice.com");
        if (!File.Exists(executable))
        {
            Console.Error.WriteLine("LibreOffice console executable not found: " + executable);
            return 2;
        }
        const string prefix = "-env:UserInstallation=file://";
        string[] normalized = args.Select(value =>
        {
            if (value.StartsWith(prefix, StringComparison.Ordinal) &&
                value.Length > prefix.Length + 2 &&
                Char.IsLetter(value[prefix.Length]) &&
                value[prefix.Length + 1] == ':' &&
                value[prefix.Length + 2] == '\\')
            {
                return prefix + "/" + value.Substring(prefix.Length).Replace('\\', '/');
            }
            return value;
        }).ToArray();
        var info = new ProcessStartInfo
        {
            FileName = executable,
            Arguments = string.Join(" ", normalized.Select(Quote)),
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            CreateNoWindow = true,
            WindowStyle = ProcessWindowStyle.Hidden,
        };
        using (var process = Process.Start(info))
        {
            string stdout = process.StandardOutput.ReadToEnd();
            string stderr = process.StandardError.ReadToEnd();
            process.WaitForExit();
            Console.Out.Write(stdout);
            Console.Error.Write(stderr);
            return process.ExitCode;
        }
    }
}
