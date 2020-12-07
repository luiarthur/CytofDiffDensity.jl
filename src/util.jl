function parse_env_var(z)
  pos = findall(r"\$\(\w+\)", z)
  if length(pos) == 0
    return z
  else
    substr = z[pos[1]]
    return ENV[substr[3:end-1]]
  end
end

function load_and_sanitize_yaml(path)
  d = YAML.load_file(path)
  return Dict((parse_env_var(k), parse_env_var(v)) for (k, v) in d)
end
